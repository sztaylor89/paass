//root script to produce 135 tof with cuts

{
  //cuts to make
  TCut vtype_cut="vtype==1";//0-small,1-medium
  TCut bsnrr_cut="bsnrr>15";
  TCut bsnrl_cut="bsnrl>15";
  TCut vsize_cut="vsize==1";
  TCut bsize_cut="bsize==1";
  TCut bqdcr_cut="bqdcr>250";
  TCut bqdcl_cut="bqdcl>250";
  TCut move_cut="move==0";
  TCut vid_cut="vid>-1";
  TCut bid_cut="bid>-1";
  TCut HPGE_cut="HPGE<99999";
  TCut cyc_cut="cyc<99999";
  TCut bcyc_cut="bcyc<99999";

  //not used: qdc, tof, snrr/l, pos, tdiff, ben, beam

  //read in bananas
  TFile *fcut = new TFile("testcut.root", "READ");
  TCutG *testcut = (TCutG*)fcut->Get("testcut");


  if(!testcut) cout << "bad\n";
  

  //hand created bananas
  //medium vandle
  //TCutG *testcut = new TCutG("testcut",11);
  testcut->SetPoint(0,90,25000);
  testcut->SetPoint(1,93,13300);
  testcut->SetPoint(2,132,8020);
  testcut->SetPoint(3,250,8020);
  testcut->SetPoint(4,250,1500);
  testcut->SetPoint(5,35,1500);
  testcut->SetPoint(6,35,5500);
  testcut->SetPoint(7,35,12000);
  testcut->SetPoint(8,35,25000);
  testcut->SetPoint(9,88,25000);
  testcut->SetPoint(10,89,25000);
  testcut->SetPoint(11,90,25000);


  //create canvas and pads
  TCanvas* c1 = new TCanvas("c1", "Tof w/ Cuts", 2000, 2000);
  gStyle->SetOptStat(111111);
  c1->Divide(2,2);
  
  
          
  //plot with cuts from above in top left
  c1->cd(1);
  V->Draw("qdc:tof>>htof(400, 0,200,5000,0,25000)",vtype_cut && bsnrr_cut && bsnrl_cut && vsize_cut && bsize_cut && bqdcr_cut && bqdcl_cut && move_cut && vid_cut && bid_cut && HPGE_cut && cyc_cut && bcyc_cut,"COLZ");

  testcut->Draw("same");

      

  //Plot banana cut in top right
  c1->cd(2);
  V->Draw("qdc:tof>>htofproj(400,0,200,5000,0,25000)",vtype_cut && bsnrr_cut && bsnrl_cut && vsize_cut && bsize_cut && bqdcr_cut && bqdcl_cut && move_cut && vid_cut && bid_cut && HPGE_cut && cyc_cut && bcyc_cut && "testcut","COLZ");

  //to use for cutting
  V->Draw("qdc:tof>>hcut(400,0,200,5000,0,25000)",vtype_cut && bsnrr_cut && bsnrl_cut && vsize_cut && bsize_cut && bqdcr_cut && bqdcl_cut && move_cut && vid_cut && bid_cut && HPGE_cut && cyc_cut && bcyc_cut && "testcut","COLZ");

  htofproj->ProjectionX()->Draw();



  //Plot shifted banana in bottom left
  //need to reset points for cut since the shift also shifts the cut
  testcut->SetPoint(0,305,25000);
  testcut->SetPoint(1,308,13300);
  testcut->SetPoint(2,347,8020);
  testcut->SetPoint(3,465,8020);
  testcut->SetPoint(4,465,1500);
  testcut->SetPoint(5,250,1500);
  testcut->SetPoint(6,250,5500);
  testcut->SetPoint(7,250,12000);
  testcut->SetPoint(8,250,25000);
  testcut->SetPoint(9,303,25000);
  testcut->SetPoint(10,304,25000);
  testcut->SetPoint(11,305,25000);

  c1->cd(3);
  V->Draw("qdc:tof-215>>htofshift(400,0,200,5000,0,25000)",vtype_cut && bsnrr_cut && bsnrl_cut && vsize_cut && bsize_cut && bqdcr_cut && bqdcl_cut && move_cut && vid_cut && bid_cut && HPGE_cut && cyc_cut && bcyc_cut && "testcut","COLZ");

  htofshift->ProjectionX()->Draw();


  //Plot banana difference in bottom right
  c1->cd(4);
  hcut->Add(htofshift,-1);
  hcut->ProjectionX()->Draw();

  //Output Raw Tof to txt
  //run convertTH1.cpp after running this script

  //save root session and pdf (might just do this via root command line)
  //just do a save as in the GUI
  
  
}
