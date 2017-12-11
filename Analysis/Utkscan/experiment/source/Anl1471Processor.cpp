/** \file Anl1471Processor.cpp
 * \brief A class to process data from ANL1471 experiment using
 * VANDLE.
 *
 *\author S. Z. Taylor
 *\date 11/5/17
 */
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

#include "BarBuilder.hpp"
#include "DammPlotIds.hpp"
#include "DoubleBetaProcessor.hpp"
#include "DetectorDriver.hpp"
#include "GeProcessor.hpp"
//#include "GetArguments.hpp"
#include "Globals.hpp"
#include "Anl1471Processor.hpp"
#include "RawEvent.hpp"
#include "TimingMapBuilder.hpp"
#include "VandleProcessor.hpp"
#include "HighResTimingData.hpp"

double Anl1471Processor::qdc_;
double Anl1471Processor::tof;


unsigned int barType;

struct Hpge_Struct{
    double ge_energy = -9998;
    double ge_time = -9998;
    double ge_beta_tdiff = -9998;
    double ge_id = -9998;
};

struct VandleRoot{
    std::vector<Hpge_Struct> HPGE;
    double tof = -9998;
    double qdc = -9998;
    double snrl = -9998;
    double snrr = -9998;
    double pos = -9998;
    double tdiff = -9998;
    double ben = -9998;
    double bqdcl = -9998;
    double bqdcr = -9998;
    double bsnrl = -9998;
    double bsnrr = -9998;
    double cyc = -9998;
    double bcyc = -9998;
    int vid = -9998;
    int vtype = -9998;
    int bid = -9998;
    int vsize = -9998;
    int bsize = -9998;

};

struct TapeInfo{
    bool move;
    bool beam;
};

struct GammaRoot{
    double gen = -9998;
    double gtime = -9998;
    double gcyc = -9998;
    double gben = -9998;
    double gbtime = -9998;
    double gbcyc = -9998;
    int gid = -9998;
    int gbid = -9998;
    int gsize = -9998;
    int bsize = -9998;
};

TapeInfo tapeinfo;
VandleRoot vroot;
GammaRoot groot;
Hpge_Struct vangaminfo;

namespace dammIds {
    namespace experiment {
        const int DD_TRACES  = 0;
        const int D_TEST  = 1;
        const int D_BADLOCATION  = 2;
        const int D_STARTLOC  = 3;
        const int DD_DEBUGGING4  = 4;
        const int DD_MedCTOFvQDC  = 5;
        const int DD_MedVetoed  = 6;
        const int DD_SmCTOFvQDC  = 7;
        const int DD_SmVetoed  = 8;
        const int DD_DEBUGGING9  = 9;
        const int D_tape = 10;
        const int D_beam = 11;
        const int DD_grow_decay = 12;
    }
}//namespace dammIds

using namespace std;
using namespace dammIds::experiment;

void Anl1471Processor::DeclarePlots(void) {
    DeclareHistogram2D(DD_TRACES, S8, SE, "traces");
    DeclareHistogram1D(D_TEST, SD, "beta gamma neutron test hist");
    DeclareHistogram1D(D_BADLOCATION, S6, "'bad' trace bar location");
    DeclareHistogram1D(D_STARTLOC, SB, "Detector Referenced as Start for vandle");
    DeclareHistogram1D(DD_DEBUGGING4, S7, "Beta Multiplicity");
    DeclareHistogram2D(DD_MedCTOFvQDC, SC, SD, "ANL-medium-<E>-vs-CorTof");
    DeclareHistogram2D(DD_MedVetoed, SC, SD, "ANL-medium-vetoed");
    DeclareHistogram2D(DD_SmCTOFvQDC, SC, SD, "ANL-small-<E>-vs-CorTof");
    DeclareHistogram2D(DD_SmVetoed, SC, SD, "ANL-small-vetoed");
    DeclareHistogram2D(DD_DEBUGGING9, SD, S6, "BSNRLvsBQDCL");
    DeclareHistogram1D(D_tape, S1, "tape move");
    DeclareHistogram1D(D_beam, S1, "beam on/off");
    DeclareHistogram2D(DD_grow_decay, SC, SA, "Grow/Decay");
}//end declare plots



Anl1471Processor::Anl1471Processor() : EventProcessor(OFFSET, RANGE, "Anl1471PRocessor") {
    associatedTypes.insert("vandle");
    associatedTypes.insert("beta");
    associatedTypes.insert("ge");

#ifdef useroot
    stringstream name;
    name << Globals::get()->GetOutputPath()
         << Globals::get()->GetOutputFileName() << ".root";
    rootfile_ = new TFile(name.str().c_str(), "RECREATE");

    roottree1_ = new TTree("V","");
    roottree2_ = new TTree("G","");

    roottree1_->Branch("vandle", &vroot, "HPGE/D:tof/D:qdc/D:snrl/D:snrr/D:pos/D:tdiff/D:ben/D:bqdcl/D:bqdcr/D:bsnrl/D:bsnrr/D:cyc/D"
            ":bcyc/D:vid/I:vtype/I:bid/I:vsize/I:bsize/I");
    roottree1_->Branch("tape", &tapeinfo,"move/b:beam/b");

    roottree2_->Branch("gamma", &groot,"gen/D:gtime/D:gcyc/D:gben/D:gbtime/D:gbcyc/D:gid/I:gbid/I:gsize/I:bsize/I");
    roottree2_->Branch("tape", &tapeinfo,"move/b:beam/b");

    QDCvsCORTOF_Medium = new TH2D("MED-QDC vs CorTof","",1100,-100,1000,25000,0,25000);
    BARvsQDC_Medium = new TH2D("MED-Bar vs QDC","",25100,-100,25000,25,0,25);
    BARvsTDIFF_Medium = new TH2D("MED-Bar vs TDIFF","",1500,-750,750,25,0,25);
    BARvsCORTOF_Medium = new TH2D("MED-Bar vs CorTof","",1100,-100,1000,25,0,25);
    QDCvsCORTOF_Small = new TH2D("SM-QDC vs CorTof","",1100,-100,1000,25000,0,25000);
    BARvsQDC_Small = new TH2D("SM-Bar vs QDC","",25100,-100,25000,25,0,25);
    BARvsTDIFF_Small = new TH2D("SM-Bar vs TDIFF","",1500,-750,750,25,0,25);
    BARvsCORTOF_Small = new TH2D("SM-Bar vs CorTof","",1100,-100,1000,25,0,25);
    GAMMA_SINGLES = new TH1D("Gamma Singles","",3100,-100,3000);
    BETA_GATED_GAMMA = new TH1D("Betad Gated Gamma","",3100,-100,3000);
    Vsize = new TH1D("Vsize","",40,0,40);
    Bsize = new TH1D("Bsize","",40,0,40);
    Gsize =new TH1D("Gsize","",40,0,40);
    BETA = new TH2D("BETA-QDCvsSNR","",8192,0,8192,64,0,64);
    GammaGrowDecay = new TH2D("Beta gated Gamma Grow/Decay","",3000,0,3000,1200,0,12);
    BetaGrowDecay = new TH2D("Beta Grow/Decay","",25000,0,25000,1200,0,12);
    NeutronGrowDecay = new TH2D("Neutron Grow/Decay","",25000,0,25000,1200,0,12);
#endif
}//end event processor





//Destructor closes/writes files and class
Anl1471Processor::~Anl1471Processor() {
#ifdef useroot
    rootfile_->Write();
    rootfile_->Close();
    delete(rootfile_);
#endif
}//end destructor



//where everything is done
bool Anl1471Processor::Process(RawEvent &event) {
    if (!EventProcessor::Process(event))
        return (false);
    //double plotMult_ = 2;
    //double plotOffset_ = 1000;

    BarMap vbars, betas;
    map<unsigned int, pair<double, double> > lrtBetas;

    BarMap betaStarts_;
    //BarMap betaSingles_;
    vector < ChanEvent * > geEvts;
    vector <vector<AddBackEvent>> geAddback;

    if (event.GetSummary("vandle")->GetList().size() != 0)
        vbars = ((VandleProcessor *) DetectorDriver::get()->
                GetProcessor("VandleProcessor"))->GetBars();


    static const vector<ChanEvent *> &doubleBetaStarts =
            event.GetSummary("beta:double:start")->GetList();
    BarBuilder startBars(doubleBetaStarts);
    startBars.BuildBars();
    betaStarts_ = startBars.GetBarMap();

    //if(event.GetSummary("beta:double")->GetList().size() != 0) {
    //    betas = ((DoubleBetaProcessor *) DetectorDriver::get()->
    //            GetProcessor("DoubleBetaProcessor"))->GetBars();
    //      if (event.GetSummary("beta:double")->GetList().size() != 0) {
    //        lrtBetas = ((DoubleBetaProcessor *) DetectorDriver::get()->
    //              GetProcessor("DoubleBetaProcessor"))->GetLowResBars();
    //}
    //}

    //static const vector<ChanEvent*> &doubleBetaSingles =
    //	event.GetSummary("beta:double:singles")->GetList();
    //BarBuilder gestartBars(doubleBetaSingles);
    //gestartBars.BuildBars();
    //betaSingles_=gestartBars.GetBarMap();


    if (event.GetSummary("ge")->GetList().size() != 0) {
        geEvts = ((GeProcessor *) DetectorDriver::get()->
                GetProcessor("GeProcessor"))->GetGeEvents();
        geAddback = ((GeProcessor *) DetectorDriver::get()->
                GetProcessor("GeProcessor"))->GetAddbackEvents();
    }

#ifdef useroot
    Vsize->Fill(vbars.size());
    Bsize->Fill(betaStarts_.size());
    Gsize->Fill(geEvts.size());


    if (TreeCorrelator::get()->place("TapeMove")->status()) {
        tapeinfo.move = 1;
    } else {
        tapeinfo.move = 0;
    }
    if (TreeCorrelator::get()->place("Beam")->status()) {
        tapeinfo.beam = 1;
    } else {
        tapeinfo.beam = 0;
    }
    plot(D_tape, tapeinfo.move);
    plot(D_beam, tapeinfo.beam);
#endif


    //plot(DD_DEBUGGING3, vbars.size());
    plot(DD_DEBUGGING4, betaStarts_.size());


    //Begin processing for VANDLE bars
    for (BarMap::iterator it = vbars.begin(); it != vbars.end(); it++) {
        TimingDefs::TimingIdentifier barId = (*it).first;
        BarDetector bar = (*it).second;

        if (!bar.GetHasEvent())
            continue;

        if (bar.GetType() == "small")
            barType = 0;
        else if (bar.GetType() == "medium")
            barType = 1;


        //stuff to test TDIFF spike
        /*if (barId.first == 2){
            plot(DD_DEBUGGING0,
                 bar.GetTimeDifference()*2+1000,
                 bar.GetLeftSide().GetMaximumValue());
            plot(DD_DEBUGGING1,
                 bar.GetTimeDifference()*2+1000,
                 bar.GetRightSide().GetMaximumValue());
            }
            plot(DD_DEBUGGING2,
             bar.GetTimeDifference()*2+1000, barId.first);
        */
        unsigned int barLoc = barId.first;
        const TimingCalibration cal = bar.GetCalibration();


        for (BarMap::iterator itStart = betaStarts_.begin();
             itStart != betaStarts_.end(); itStart++) {
            BarDetector beta_start = (*itStart).second;
            unsigned int startLoc = (*itStart).first.first;
            if (!beta_start.GetHasEvent())
                continue;
            double tofOffset = cal.GetTofOffset(startLoc);
            double tof = bar.GetCorTimeAve() -
                         beta_start.GetCorTimeAve() + tofOffset;


            double corTof =
                    ((VandleProcessor *) DetectorDriver::get()->
                            GetProcessor("VandleProcessor"))->
                            CorrectTOF(tof, bar.GetFlightPath(), cal.GetZ0());


            //tape move veto cut damm
            bool tapeMove = TreeCorrelator::get()->place("TapeMove")->status();
            if (tapeMove == 0) { //plot only if tape is NOT moving
                if (bar.GetType() == "medium")
                    plot(DD_MedCTOFvQDC, corTof * 2 + 1000, bar.GetQdc());

                if (bar.GetType() == "small")
                    plot(DD_SmCTOFvQDC, corTof * 2 + 1000, bar.GetQdc());
            }

            if (tapeMove == 1) { //plot only if tape is moving
                if (bar.GetType() == "medium")
                    plot(DD_MedVetoed, corTof * 2 + 1000, bar.GetQdc());

                if (bar.GetType() == "small")
                    plot(DD_SmVetoed, corTof * 2 + 1000, bar.GetQdc());
            }

            plot(DD_DEBUGGING9, beta_start.GetLeftSide().GetTraceQdc(),
                 beta_start.GetLeftSide().GetTrace().GetSignalToNoiseRatio());

            //adding cycle time stuff for vandle
            double vcyc_time, bcyc_time, cyc_time, vtime, btime;
            vcyc_time = bcyc_time = cyc_time = vtime = btime = -9999;
            if (TreeCorrelator::get()->place("Cycle")->status()) {
                cyc_time = TreeCorrelator::get()->place("Cycle")->last().time;
                cyc_time *= (Globals::get()->GetClockInSeconds() * 1.e9);//converts from clock ticks to ns

                //vandle event
                vtime = bar.GetCorTimeAve();//in ns
                vcyc_time = (vtime - cyc_time);//in ns
                vcyc_time /= 1e9;//converts from ns to s

                //beta event
                btime = beta_start.GetCorTimeAve();//in ns
                bcyc_time = (btime - cyc_time);//in ns
                bcyc_time /= 1e9;//converts from ns to s
            }

            //adding HPGE energy info to vandle tree
            if (geEvts.size() != 0) {
                for (vector<ChanEvent *>::const_iterator itHPGE = geEvts.begin();
                     itHPGE != geEvts.end(); itHPGE++){
                    double B_time, G_time;
                    vangaminfo.ge_id = (*itHPGE)->GetChanID().GetLocation();
                    G_time = (*itHPGE)->GetTimeSansCfd();//gives result in clock ticks
                    //double GG = (*itHPGE)->GetTimeSansCfd();// used as check
                    G_time *= Globals::get()->GetClockInSeconds() * 1.e9; //converts clock ticks to ns
                    B_time = beta_start.GetCorTimeAve(); //gives result in ns
                    vangaminfo.ge_beta_tdiff = G_time - B_time;
                    vangaminfo.ge_time=G_time;
                    vangaminfo.ge_energy = (*itHPGE)->GetCalibratedEnergy();
                    vroot.HPGE.emplace_back(vangaminfo);
                }
            }





#ifdef useroot
            vroot.tof = corTof;
            vroot.qdc = bar.GetQdc();
            vroot.snrl = bar.GetLeftSide().GetTrace().GetSignalToNoiseRatio();
            vroot.snrr = bar.GetRightSide().GetTrace().GetSignalToNoiseRatio();
            vroot.pos = bar.GetQdcPosition();
            vroot.tdiff = bar.GetTimeDifference();
            vroot.ben = beta_start.GetQdc();
            vroot.bqdcl = beta_start.GetLeftSide().GetTraceQdc();
            vroot.bqdcr = beta_start.GetRightSide().GetTraceQdc();
            vroot.bsnrl = beta_start.GetLeftSide().GetTrace().GetSignalToNoiseRatio();
            vroot.bsnrr = beta_start.GetRightSide().GetTrace().GetSignalToNoiseRatio();
            vroot.cyc = vcyc_time;
            vroot.bcyc = bcyc_time;
            //vroot.HPGE = vangaminfo;  May not need this, emplace may fill this above.
            vroot.vid = barLoc;
            vroot.vtype = barType;
            vroot.bid = startLoc;
            vroot.vsize = vbars.size();
            vroot.bsize = betaStarts_.size();

#endif


#ifdef useroot
            if (barType == 1) {
                QDCvsCORTOF_Medium->Fill(corTof, bar.GetQdc());
                BARvsQDC_Medium->Fill(bar.GetQdc(), barLoc);
                BARvsTDIFF_Medium->Fill(bar.GetTimeDifference(), barLoc);
                BARvsCORTOF_Medium->Fill(corTof, barLoc);
            } else if (barType == 0) {
                QDCvsCORTOF_Small->Fill(corTof, bar.GetQdc());
                BARvsQDC_Small->Fill(bar.GetQdc(), barLoc);
                BARvsTDIFF_Small->Fill(bar.GetTimeDifference(), barLoc);
                BARvsCORTOF_Small->Fill(corTof, barLoc);
            }
            BETA->Fill(vroot.bqdcl, vroot.bsnrl);
            BetaGrowDecay->Fill(beta_start.GetQdc(),bcyc_time);
            NeutronGrowDecay->Fill(bar.GetQdc(),vcyc_time);
            qdc_ = bar.GetQdc();
            //tof = tof;
            roottree1_->Fill();

            //clear vector and "zero" structure for gamma stuff
            vroot.HPGE.clear();
            vangaminfo.ge_beta_tdiff=-9999;
            vangaminfo.ge_energy=-9999;
            vangaminfo.ge_id=-9999;
            vangaminfo.ge_time=-9999;
            // bar.GetLeftSide().ZeroRootStructure(leftVandle);
            // bar.GetRightSide().ZeroRootStructure(rightVandle);
            // beta_start.GetLeftSide().ZeroRootStructure(leftBeta);
            // beta_start.GetRightSide().ZeroRootStructure(rightBeta);
            qdc_ = tof = -9999;
            //VID = BID = SNRVL = SNRVR = -9999;
            //GamEn = SNRBL = SNRBR = vandle_ = beta_ = ge_ = -9999;
#endif

            //            plot(DD_DEBUGGING1, tof * plotMult_ + plotOffset_, bar.GetQdc());

        } // for(TimingMap::iterator itStart
    } //(BarMap::iterator itBar
    //End processing for VANDLE bars





    //Stuff to fill HPGe branch

    if (geEvts.size() != 0) {
        for (vector<ChanEvent *>::const_iterator itGe = geEvts.begin();
             itGe != geEvts.end(); itGe++) {
            double ge_energy, ge_time, gb_time, grow_decay_time, gb_en, gcyc_time, gb_grow_decay_time;
            ge_energy = ge_time = gb_time = grow_decay_time = gb_en = gcyc_time = gb_grow_decay_time =-9999.0;
            int ge_id = -9999;
            int gb_startLoc = -9999;
            BarDetector gb_start;
            ge_energy = (*itGe)->GetCalibratedEnergy();
            ge_id = (*itGe)->GetChanID().GetLocation();
            ge_time = (*itGe)->GetWalkCorrectedTime();
            ge_time *= (Globals::get()->GetClockInSeconds() * 1.e9);//converts from clock ticks to ns

            if (TreeCorrelator::get()->place("Cycle")->status()) {
                gcyc_time = TreeCorrelator::get()->place("Cycle")->last().time;
                gcyc_time *= (Globals::get()->GetClockInSeconds() * 1.e9);//converts from clock ticks to ns
                grow_decay_time = (ge_time - gcyc_time);//in ns
                grow_decay_time /= 1e9;//converts from ns to s
                //cout << ge_energy << endl << grow_decay_time << endl << endl;
                //plot(DD_grow_decay, ge_energy, grow_decay_time);
            }

            if (doubleBetaStarts.size() != 0) {
                for (BarMap::iterator itGB = betaStarts_.begin();
                     itGB != betaStarts_.end(); itGB++) {
                    gb_start = (*itGB).second;
                    gb_startLoc = (*itGB).first.first;
                    gb_en = gb_start.GetQdc();
                    gb_time = gb_start.GetCorTimeAve();//in ns
                    gb_grow_decay_time = (gb_time - gcyc_time);//in ns
                    gb_grow_decay_time /= 1e9;//converts from ns to s
                }
            } else {
                gb_startLoc = -8888;
                gb_en = -8888;
                gb_time = -8888;
            }

#ifdef useroot
            groot.gen = ge_energy;
            groot.gtime = ge_time;
            groot.gcyc = grow_decay_time;
            groot.gben = gb_en;
            groot.gbtime = gb_time;
            groot.gbcyc = gb_grow_decay_time;
            groot.gid = ge_id;
            groot.gbid = gb_startLoc;
            groot.gsize = geEvts.size();
            groot.bsize = betaStarts_.size();

            roottree2_->Fill();
            GAMMA_SINGLES->Fill(ge_energy);
            //            GrowDecay->Fill(ge_energy,grow_decay_time);
            if (doubleBetaStarts.size() != 0) {
                BETA_GATED_GAMMA->Fill(ge_energy);
                GammaGrowDecay->Fill(ge_energy,grow_decay_time);
                plot(DD_grow_decay, ge_energy, grow_decay_time);
            }
#endif
        }
    }


    EndProcess();
    return(true);
}