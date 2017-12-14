/** \file Anl1471Processor.cpp
 * \brief A class to process data from ANL1471 experiment using VANDLE.
 * \brief This class creates 2 root trees, one for neutron based events and the other for gamma based events.
 *\author S. Z. Taylor
 *\date 11/5/17
 */
#include <fstream>
#include <iostream>
#include <cmath>

#include "BarBuilder.hpp"
#include "DammPlotIds.hpp"
#include "DoubleBetaProcessor.hpp"
#include "DetectorDriver.hpp"
#include "GeProcessor.hpp"
#include "Anl1471Processor.hpp"
#include "VandleProcessor.hpp"

//double Anl1471Processor::tof;


//neutron related variables
std::vector<double> tof;//time of flight
std::vector<double> qdc;//vandle vbarenergy/qdc
std::vector<double> snrl;//vandle signal to noise ratio left
std::vector<double> snrr;//vandle signal to noise ratio right
std::vector<double> pos;//position in bar
std::vector<double> tdiff;//time difference between vandle bar ends
std::vector<double> ben;//beta bar energy/qdc
std::vector<double> bqdcl;//beta qdc left
std::vector<double> bqdcr;//beta qdc right
std::vector<double> bsnrl;// beta signal to noise ratio left
std::vector<double> bsnrr;//beta signal to noise ratio right
std::vector<double> cyc;//vandle event time in reference to start of cycle
std::vector<double> bcyc;//beta event time in reference to  start of cycle
std::vector<double> gen;//gamma energy
std::vector<double> Vtime;//vandle event time
std::vector<double> Gtime;//gamma event time
std::vector<double> Btime;//beta event time
std::vector<double> BGtime;//time difference between beta and gamma events
std::vector<int> vid;//vandle bar id
std::vector<int> vtype;//vandle bar type, 0=small, 1=medium
std::vector<int> bid;//beta bar id
std::vector<int> gid;//hpge clover leaf id
std::vector<int> vsize;//size of vandle event
std::vector<int> gsize;//size of gamma event
std::vector<int> bsize;//size of beta event

//tape related variables
std::vector<bool> move;//moving tape collector moving/stationary   0=stopped, 1=moving
std::vector<bool> beam;//beam on/off    0=off, 1=on

//gamma singles variables, keeping separate vectors to avoid confusion from similar ones in neutron variables
std::vector<double> g_en;//gamma energy
std::vector<double> g_time;//gamma event time
std::vector<double> g_cyc;//gamma event time in reference to start of cycle
std::vector<double> g_ben;//beta energy
std::vector<double> g_btime;//beta event time
std::vector<double> g_bcyc;//beta event time in reference to start of cycle
std::vector<int> g_id;//hpge clover leaf id
std::vector<int> g_bid;//beta bar id
std::vector<int> g_size;//size of gamma event
std::vector<int> g_bsize;//size of beta event

//assigning damm histogram id's,    Starts at 6050
namespace dammIds {
    namespace experiment {
     /*   const int DD_TRACES  = 0;
        const int D_TEST  = 1;
        const int D_BADLOCATION  = 2;
        const int D_STARTLOC  = 3;
        */
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

//declaring damm plots
void Anl1471Processor::DeclarePlots() {
/*    DeclareHistogram2D(DD_TRACES, S8, SE, "traces");
    DeclareHistogram1D(D_TEST, SD, "beta gamma neutron test hist");
    DeclareHistogram1D(D_BADLOCATION, S6, "'bad' trace bar location");
    DeclareHistogram1D(D_STARTLOC, SB, "Detector Referenced as Start for vandle");
    DeclareHistogram1D(DD_DEBUGGING4, S7, "Beta Multiplicity");
*/
    DeclareHistogram2D(DD_MedCTOFvQDC, SC, SD, "ANL-medium-<E>-vs-CorTof");
    DeclareHistogram2D(DD_MedVetoed, SC, SD, "ANL-medium-vetoed");
    DeclareHistogram2D(DD_SmCTOFvQDC, SC, SD, "ANL-small-<E>-vs-CorTof");
    DeclareHistogram2D(DD_SmVetoed, SC, SD, "ANL-small-vetoed");
    DeclareHistogram2D(DD_DEBUGGING9, SD, S6, "BSNRLvsBQDCL");
    DeclareHistogram1D(D_tape, S1, "tape move");
    DeclareHistogram1D(D_beam, S1, "beam on/off");
    DeclareHistogram2D(DD_grow_decay, SC, SA, "Grow/Decay");
}//end declare plots


//gathering detector types
Anl1471Processor::Anl1471Processor() : EventProcessor(OFFSET, RANGE, "Anl1471PRocessor") {
    associatedTypes.insert("vandle");
    associatedTypes.insert("beta");
    associatedTypes.insert("ge");

    //automatically creating root file with same name as damm output
#ifdef useroot
    stringstream name;
    name << Globals::get()->GetOutputPath()
         << Globals::get()->GetOutputFileName() << ".root";
    rootfile_ = new TFile(name.str().c_str(), "RECREATE");

    //creating 2 root trees, one based of of neutron events(b-n and b-n-y), the other off of gamma events(y singles and b-y)
    roottree1_ = new TTree("V","Tree with neutron related events");
    roottree2_ = new TTree("G","Tree with gamma related events");

    //defining branches for each vector variable
    //neutron based
    roottree1_->Branch("tof",&tof);
    roottree1_->Branch("qdc",&qdc);
    roottree1_->Branch("snrl",&snrl);
    roottree1_->Branch("snrr",&snrr);
    roottree1_->Branch("pos",&pos);
    roottree1_->Branch("tdiff",&tdiff);
    roottree1_->Branch("ben",&ben);
    roottree1_->Branch("bqdcl",&bqdcl);
    roottree1_->Branch("bqdcr",&bqdcr);
    roottree1_->Branch("bsnrl",&bsnrl);
    roottree1_->Branch("bsnrr",&bsnrr);
    roottree1_->Branch("cyc",&cyc);
    roottree1_->Branch("bcyc",&bcyc);
    roottree1_->Branch("gen",&gen);
    roottree1_->Branch("Vtime",&Vtime);
    roottree1_->Branch("Gtime",&Gtime);
    roottree1_->Branch("Btime",&Btime);
    roottree1_->Branch("BGtime",&BGtime);
    roottree1_->Branch("vid",&vid);
    roottree1_->Branch("vtype",&vtype);
    roottree1_->Branch("bid",&bid);
    roottree1_->Branch("gid",&gid);
    roottree1_->Branch("vsize",&vsize);
    roottree1_->Branch("gsize",&gsize);
    roottree1_->Branch("bsize",&bsize);
    roottree1_->Branch("move",&move);
    roottree1_->Branch("beam",&beam);

    //gamma based
    roottree2_->Branch("g_en",&g_en);
    roottree2_->Branch("g_time",&g_time);
    roottree2_->Branch("g_cyc",&g_cyc);
    roottree2_->Branch("g_ben",&g_ben);
    roottree2_->Branch("g_btime",&g_btime);
    roottree2_->Branch("g_bcyc",&g_bcyc);
    roottree2_->Branch("g_id",&g_id);
    roottree2_->Branch("g_bid",&g_bid);
    roottree2_->Branch("g_size",&g_size);
    roottree2_->Branch("g_bsize",&g_bsize);
    roottree2_->Branch("move",&move);
    roottree2_->Branch("beam",&beam);

    //pre-making some frequently used histograms in root
    QDCvsCORTOF_Medium = new TH2D("MED-QDC vs CorTof","",1100,-100,1000,25000,0,25000);
    BARvsQDC_Medium = new TH2D("MED-Bar vs QDC","",25100,-100,25000,25,0,25);
    BARvsTDIFF_Medium = new TH2D("MED-Bar vs TDIFF","",200,-25,25,25,0,25);
    BARvsCORTOF_Medium = new TH2D("MED-Bar vs CorTof","",1100,-100,1000,25,0,25);
    QDCvsCORTOF_Small = new TH2D("SM-QDC vs CorTof","",1100,-100,1000,25000,0,25000);
    BARvsQDC_Small = new TH2D("SM-Bar vs QDC","",25100,-100,25000,25,0,25);
    BARvsTDIFF_Small = new TH2D("SM-Bar vs TDIFF","",200,-25,25,25,0,25);
    BARvsCORTOF_Small = new TH2D("SM-Bar vs CorTof","",1100,-100,1000,25,0,25);
    GAMMA_SINGLES = new TH1D("Gamma Singles","",3100,-100,3000);
    BETA_GATED_GAMMA = new TH1D("Beta Gated Gamma","",3100,-100,3000);
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

//process where all calculations and filling occurs
bool Anl1471Processor::Process(RawEvent &event) {
    if (!EventProcessor::Process(event))
        return (false);

    //creating and getting detectors
    BarMap vbars, betaStarts;
    vector < ChanEvent * > geEvts;
    vector <vector<AddBackEvent>> geAddback;

    if (event.GetSummary("vandle")->GetList().empty()) {
        vbars = ((VandleProcessor *) DetectorDriver::get()->
                GetProcessor("VandleProcessor"))->GetBars();
    }

    static const vector<ChanEvent *> &doubleBetaStarts =
            event.GetSummary("beta:double:start")->GetList();
    BarBuilder startBars(doubleBetaStarts);
    startBars.BuildBars();
    betaStarts = startBars.GetBarMap();


    if (event.GetSummary("ge")->GetList().empty()) {
        geEvts = ((GeProcessor *) DetectorDriver::get()->
                GetProcessor("GeProcessor"))->GetGeEvents();
        geAddback = ((GeProcessor *) DetectorDriver::get()->
                GetProcessor("GeProcessor"))->GetAddbackEvents();
    }

#ifdef useroot
    //fill pre-made size histograms
    Vsize->Fill(vbars.size());
    Bsize->Fill(betaStarts.size());
    Gsize->Fill(geEvts.size());

//filling tape related vectors
//clearing first to make sure it is empty
    move.clear();
    beam.clear();
    //varibles used to fill vectors and damm plots
    int MOVE, BEAM;

    if (TreeCorrelator::get()->place("TapeMove")->status()) {
        move.emplace_back(MOVE);
    } else {
        move.emplace_back(MOVE);
    }

    if (TreeCorrelator::get()->place("Beam")->status()) {
        beam.emplace_back(BEAM);
    } else {
        beam.emplace_back(BEAM);
    }
    plot(D_tape, MOVE);
    plot(D_beam, BEAM);
#endif


    plot(DD_DEBUGGING4, betaStarts.size());


    //Begin processing for VANDLE bars
    for (BarMap::iterator it = vbars.begin(); it != vbars.end(); it++) {
        TimingDefs::TimingIdentifier barId = (*it).first;
        BarDetector bar = (*it).second;
        int barType = -9999;

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


        for (BarMap::iterator itStart = betaStarts.begin();
             itStart != betaStarts.end(); itStart++) {
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
            double HPGE_energy = -9999.0;
            double BG_TDIFF = -9999.0;
            int gamma_id=-9999;
            if (geEvts.size() != 0) {
	      for (vector<ChanEvent *>::const_iterator itHPGE = geEvts.begin();
		   itHPGE != geEvts.end(); itHPGE++){
                double B_time, G_time;
                gamma_id = (*itHPGE)->GetChanID().GetLocation();
                G_time = (*itHPGE)->GetTimeSansCfd();//gives result in clock ticks
                //double GG = (*itHPGE)->GetTimeSansCfd();// used as check
                G_time *= Globals::get()->GetClockInSeconds() * 1.e9; //converts clock ticks to ns
                B_time = beta_start.GetCorTimeAve(); //gives result in ns
                BG_TDIFF = G_time - B_time;
                if (BG_TDIFF > 0){
                    HPGE_energy = (*itHPGE)->GetCalibratedEnergy();
                    plot (D_TEST, HPGE_energy);
                }else {
                    HPGE_energy = -7777.0;
                    gamma_id = -7777;
                }
          }
	    }else{ 
                HPGE_energy = -8888.0;
                gamma_id=-8888;
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
            vroot.HPGE = HPGE_energy;
            vroot.BGtime = BG_TDIFF;
            vroot.vid = barLoc;
            vroot.vtype = barType;
            vroot.bid = startLoc;
            vroot.gid = gamma_id;
            vroot.vsize = vbars.size();
            vroot.bsize = betaStarts.size();

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
            roottree1_->Fill();

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
                for (BarMap::iterator itGB = betaStarts.begin();
                     itGB != betaStarts.end(); itGB++) {
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
            groot.bsize = betaStarts.size();

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

