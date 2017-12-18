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
std::vector<int> v_move;//moving tape collector moving/stationary   0=stopped, 1=moving
std::vector<int> v_beam;//beam on/off    0=off, 1=on
std::vector<int> g_move;//move for gamma tree
std::vector<int> g_beam;//beam for gamma tree

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

//assigning damm histogram id's,    Starts at damm id 6050
namespace dammIds {
    namespace experiment {
        const int DD_MedCTOFvQDC  = 0;
        const int DD_MedVetoed  = 1;
        const int DD_SmCTOFvQDC  = 2;
        const int DD_SmVetoed  = 3;
        const int DD_SNRvQDC  = 4;
        const int D_tape = 5;
        const int D_beam = 6;
        const int DD_grow_decay = 7;
    }
}//namespace dammIds

using namespace std;
using namespace dammIds::experiment;

//declaring damm plots
void Anl1471Processor::DeclarePlots() {
    DeclareHistogram2D(DD_MedCTOFvQDC, SC, SD, "ANL-medium-<E>-vs-CorTof");
    DeclareHistogram2D(DD_MedVetoed, SC, SD, "ANL-medium-vetoed");
    DeclareHistogram2D(DD_SmCTOFvQDC, SC, SD, "ANL-small-<E>-vs-CorTof");
    DeclareHistogram2D(DD_SmVetoed, SC, SD, "ANL-small-vetoed");
    DeclareHistogram2D(DD_SNRvQDC, SD, S6, "BSNRLvsBQDCL");
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
    roottree1_->Branch("move",&v_move);
    roottree1_->Branch("beam",&v_beam);

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
    roottree2_->Branch("move",&g_move);
    roottree2_->Branch("beam",&g_beam);

    //pre-making some frequently used histograms in root
    QDCvsCORTOF_Medium = new TH2D("MED-QDC vs CorTof","",1100,-100.0,1000,25000,0,25000);
    BARvsQDC_Medium = new TH2D("MED-Bar vs QDC","",25100,-100.0,25000,25,0,25);
    BARvsTDIFF_Medium = new TH2D("MED-Bar vs TDIFF","",200,-25.0,25,25,0,25);
    BARvsCORTOF_Medium = new TH2D("MED-Bar vs CorTof","",1100,-100.0,1000,25,0,25);
    QDCvsCORTOF_Small = new TH2D("SM-QDC vs CorTof","",1100,-100.0,1000,25000,0,25000);
    BARvsQDC_Small = new TH2D("SM-Bar vs QDC","",25100,-100.0,25000,25,0,25);
    BARvsTDIFF_Small = new TH2D("SM-Bar vs TDIFF","",200,-25.0,25,25,0,25);
    BARvsCORTOF_Small = new TH2D("SM-Bar vs CorTof","",1100,-100.0,1000,25,0,25);
    GAMMA_SINGLES = new TH1D("Gamma Singles","",3100,-100.0,3000);
    BETA_GATED_GAMMA = new TH1D("Beta Gated Gamma","",3100,-100.0,3000);
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
        vbars = ((VandleProcessor *) DetectorDriver::get()->GetProcessor("VandleProcessor"))->GetBars();
    }

    static const vector<ChanEvent *> &doubleBetaStarts = event.GetSummary("beta:double:start")->GetList();
    BarBuilder startBars(doubleBetaStarts);
    startBars.BuildBars();
    betaStarts = startBars.GetBarMap();


    if (event.GetSummary("ge")->GetList().empty()) {
        geEvts = ((GeProcessor *) DetectorDriver::get()->GetProcessor("GeProcessor"))->GetGeEvents();
        geAddback = ((GeProcessor *) DetectorDriver::get()->GetProcessor("GeProcessor"))->GetAddbackEvents();
    }

#ifdef useroot
    //fill pre-made size histograms
    Vsize->Fill(vbars.size());
    Bsize->Fill(betaStarts.size());
    Gsize->Fill(geEvts.size());

//filling tape related vectors
//clearing first to make sure it is empty
    v_move.clear();
    v_beam.clear();
    g_move.clear();
    g_beam.clear();

    //variables used to fill vectors and damm plots
    int MOVE, BEAM;

    if (TreeCorrelator::get()->place("TapeMove")->status()) {
        MOVE = 1;
        v_move.emplace_back(MOVE);
        g_move.emplace_back(MOVE);
    } else {
        MOVE = 0;
        v_move.emplace_back(MOVE);
        g_move.emplace_back(MOVE);
    }

    if (TreeCorrelator::get()->place("Beam")->status()) {
        BEAM = 1;
        v_beam.emplace_back(BEAM);
        g_beam.emplace_back(BEAM);
    } else {
        BEAM = 0;
        v_beam.emplace_back(BEAM);
        g_beam.emplace_back(BEAM);
    }
#endif

    //plot tape related damm histograms
    plot(D_tape, MOVE);
    plot(D_beam, BEAM);


    //Begin processing for VANDLE bars
    //loop over vandle bar events
    for (BarMap::iterator it = vbars.begin(); it != vbars.end(); it++) {

        //getting info from iterator
        TimingDefs::TimingIdentifier barId = (*it).first;
        BarDetector bar = (*it).second;
        int barType = -9999;
        unsigned int barLoc = barId.first;
        const TimingCalibration cal = bar.GetCalibration();

        //check in case code didn't have event
        if (!bar.GetHasEvent())
            continue;

        //set bar type
        if (bar.GetType() == "small")
            barType = 0;
        else if (bar.GetType() == "medium")
            barType = 1;

        //loop over beta events
        for (BarMap::iterator itStart = betaStarts.begin(); itStart != betaStarts.end(); itStart++) {

            //clearing vectors.  Since I fill my root tree in the beta loop I will clear vectors here
            tof.clear();
            qdc.clear();
            snrl.clear();
            snrr.clear();
            pos.clear();
            tdiff.clear();
            ben.clear();
            bqdcl.clear();
            bqdcr.clear();
            bsnrl.clear();
            bsnrr.clear();
            cyc.clear();
            bcyc.clear();
            gen.clear();
            Vtime.clear();
            Gtime.clear();
            Btime.clear();
            BGtime.clear();
            vid.clear();
            vtype.clear();
            bid.clear();
            gid.clear();
            vsize.clear();
            gsize.clear();
            bsize.clear();

            BarDetector beta_start = (*itStart).second;
            unsigned startLoc = (*itStart).first.first;

            //check for event
            if (!beta_start.GetHasEvent())
                continue;

            double tofOffset = cal.GetTofOffset(startLoc);
            double TOF = bar.GetCorTimeAve() - beta_start.GetCorTimeAve() + tofOffset;
            double corTof =
                    ((VandleProcessor *) DetectorDriver::get()->GetProcessor("VandleProcessor"))->CorrectTOF(TOF, bar
                            .GetFlightPath(), cal.GetZ0());

            //tape move veto cut for damm
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

            plot(DD_SNRvQDC, beta_start.GetLeftSide().GetTraceQdc(), beta_start.GetLeftSide().GetTrace().GetSignalToNoiseRatio());

            double vcyc_time, bcyc_time, btime, vtime;
            vcyc_time = bcyc_time = btime = vtime = -9999;
            //cycle time stuff for vandle
            if (TreeCorrelator::get()->place("Cycle")->status()) {
                double cyc_time = TreeCorrelator::get()->place("Cycle")->last().time;
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
            unsigned int gamma_id = 9999;
            //get hpge event info if it exists
            if (geEvts.size() != 0) {
                //loop over hpge events
                for (vector<ChanEvent *>::const_iterator itHPGE = geEvts.begin(); itHPGE != geEvts.end(); itHPGE++){
 //                   gamma_id = (*itHPGE)->GetChanID().GetLocation(); //ERROR
                    double G_time = (*itHPGE)->GetTimeSansCfd();//gives result in clock ticks
                    G_time *= Globals::get()->GetClockInSeconds() * 1.e9; //converts clock ticks to ns
                    double B_time = beta_start.GetCorTimeAve(); //gives result in ns
                    BG_TDIFF = G_time - B_time;
                    HPGE_energy = (*itHPGE)->GetCalibratedEnergy();

                    //fill hpge vectors
                    gen.emplace_back(HPGE_energy);
                    Gtime.emplace_back(G_time);
                    BGtime.emplace_back(BG_TDIFF);
                    gid.emplace_back(gamma_id);
                    gsize.emplace_back(geEvts.size());

                }
            }else{//fill hpge vectors with default value to cut out later in root
                gen.emplace_back(-8888);
                Gtime.emplace_back(-8888);
                BGtime.emplace_back(-8888);
                gid.emplace_back(-8888);
                gsize.emplace_back(-8888);
            }

//fill vandle and beta vectors
#ifdef useroot
            tof.emplace_back(corTof);
            qdc.emplace_back(bar.GetQdc());
            snrl.emplace_back(bar.GetLeftSide().GetTrace().GetSignalToNoiseRatio());
            snrr.emplace_back(bar.GetRightSide().GetTrace().GetSignalToNoiseRatio());
            pos.emplace_back(bar.GetQdcPosition());
            tdiff.emplace_back(bar.GetTimeDifference());
            ben.emplace_back(beta_start.GetQdc());
            bqdcl.emplace_back(beta_start.GetLeftSide().GetTraceQdc());
            bqdcr.emplace_back(beta_start.GetRightSide().GetTraceQdc());
            bsnrl.emplace_back(beta_start.GetLeftSide().GetTrace().GetSignalToNoiseRatio());
            bsnrr.emplace_back(beta_start.GetRightSide().GetTrace().GetSignalToNoiseRatio());
            cyc.emplace_back(vcyc_time);
            bcyc.emplace_back(bcyc_time);
            Vtime.emplace_back(vtime);
            Btime.emplace_back(btime);
            vid.emplace_back(barLoc);
            vtype.emplace_back(barType);
            bid.emplace_back(startLoc);
            vsize.emplace_back(vbars.size());
            bsize.emplace_back(betaStarts.size());

            //fill pre-made vandle related root histograms
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
            BETA->Fill(beta_start.GetLeftSide().GetTraceQdc(), beta_start.GetLeftSide().GetTrace().GetSignalToNoiseRatio());//bqdcl vs bsnrl
            BetaGrowDecay->Fill(beta_start.GetQdc(),bcyc_time);
            NeutronGrowDecay->Fill(bar.GetQdc(),vcyc_time);

            //fill vandle root tree
            roottree1_->Fill();
#endif
        } // for(TimingMap::iterator itStart
    } //(BarMap::iterator itBar
    //End processing for VANDLE bars





    //Stuff to fill HPGe branch

    if (geEvts.size() != 0) {
        for (vector<ChanEvent *>::const_iterator itGe = geEvts.begin(); itGe != geEvts.end(); itGe++) {

            //clear vectors here since this for loop fills the root tree
            g_en.clear();
            g_time.clear();
            g_cyc.clear();
            g_ben.clear();
            g_btime.clear();
            g_bcyc.clear();
            g_id.clear();
            g_bid.clear();
            g_size.clear();
            g_bsize.clear();

            //variables used for calculations and filling vectors
            double gb_time, grow_decay_time, gb_en, gcyc_time, gb_grow_decay_time;
            gb_time = grow_decay_time = gb_en = gcyc_time = gb_grow_decay_time =-9999.0;
            int gb_startLoc = -9999;

            //creating beta bar
            BarDetector gb_start;
            double ge_energy = (*itGe)->GetCalibratedEnergy();
//            unsigned int ge_id = (*itGe)->GetChanID().GetLocation(); //ERROR
            unsigned int ge_id = 0;
            double ge_time = (*itGe)->GetWalkCorrectedTime();
            ge_time *= (Globals::get()->GetClockInSeconds() * 1.e9);//converts from clock ticks to ns

            if (TreeCorrelator::get()->place("Cycle")->status()) {
                gcyc_time = TreeCorrelator::get()->place("Cycle")->last().time;
                gcyc_time *= (Globals::get()->GetClockInSeconds() * 1.e9);//converts from clock ticks to ns
                grow_decay_time = (ge_time - gcyc_time);//in ns
                grow_decay_time /= 1e9;//converts from ns to s
            }

            if (doubleBetaStarts.size() != 0) {
                for (BarMap::iterator itGB = betaStarts.begin(); itGB != betaStarts.end(); itGB++) {
                    gb_start = (*itGB).second;
                    gb_startLoc = (*itGB).first.first;
                    gb_en = gb_start.GetQdc();
                    gb_time = gb_start.GetCorTimeAve();//in ns
                    gb_grow_decay_time = (gb_time - gcyc_time);//in ns
                    gb_grow_decay_time /= 1e9;//converts from ns to s

                    //fill beta vectors
                    g_ben.emplace_back(gb_en);
                    g_btime.emplace_back(gb_time);
                    g_bcyc.emplace_back(gb_grow_decay_time);
                    g_bid.emplace_back(gb_startLoc);
                    g_bsize.emplace_back(betaStarts.size());
                }
            } else {//fill beta vectors with default value to cut out later in root
                g_ben.emplace_back(-8888);
                g_btime.emplace_back(-8888);
                g_bcyc.emplace_back(-8888);
                g_bid.emplace_back(-8888);
                g_bsize.emplace_back(-8888);
            }

#ifdef useroot
            //fill gamma vectors
            g_en.emplace_back(ge_energy);
            g_time.emplace_back(ge_time);
            g_cyc.emplace_back(grow_decay_time);
            g_id.emplace_back(ge_id);
            g_size.emplace_back(geEvts.size());

            //fill pre-made gamma root histograms
            GAMMA_SINGLES->Fill(ge_energy);
            if (doubleBetaStarts.size() != 0) {
                BETA_GATED_GAMMA->Fill(ge_energy);
                GammaGrowDecay->Fill(ge_energy,grow_decay_time);
                plot(DD_grow_decay, ge_energy, grow_decay_time);
            }

            //fill gamma root tree
            roottree2_->Fill();
#endif
        }//for itGe
    }//if geEvts

    EndProcess();
    return(true);
}

