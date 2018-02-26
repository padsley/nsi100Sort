#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>

TFile *fin, *frawin;
TTree *trin, *traw;

void OpenFile()
{
  fin = TFile::Open("output.root");
  trin = (TTree*)fin->FindObjectAny("DATA");
  //trin->Print();
}		   

void OpenRawFile(int RunNumber)
{
  char buffer[256];
  sprintf(buffer,"../sorted/R%d_0_old.root",RunNumber);
  frawin = TFile::Open(buffer);
  traw = (TTree*)frawin->FindObjectAny("EGTree");
}

void CheckRun(int RunNumber)
{
  OpenFile();
  OpenRawFile(RunNumber);
  TCanvas *RunCanvas = new TCanvas("RunCanvas","RunCanvas",1800,1600);

  RunCanvas->Divide(3,2);

  RunCanvas->cd(1);
  char buffer[256];

  //sprintf(buffer,"run==%d",RunNumber);

  traw->Draw("tdcList:adcList>>hADCTDC(192,0,192,112,0,112)","","col");


  RunCanvas->cd(2);
  int nentries = traw->GetEntries();
  sprintf(buffer,"SPpos:Entry$>>hPosEntry(%d,0,%d,512,0,4096)",nentries/100,nentries);
  traw->Draw(buffer,"","col");
  hPosEntry->SetStats(0);

  RunCanvas->cd(3);
  sprintf(buffer,"run==%d",RunNumber);
  trin->Draw("SiliconTDCChannel:TimeDifference>>hTDCPerChannel(4000,-1000,3000,96,0,96)",buffer,"col");
  //  hTDCPerChannel->GetXaxis()->SetRangeUser(-1000,1500);
  TPad *paddy3 = (TPad*)RunCanvas->FindObject("RunCanvas_3");
  paddy3->SetLogz();
  

  RunCanvas->cd(4);
  traw->Draw("tdcList:tdcData>>h2DTDCValueVsTDCChannel(2500,0,2500,112,0,112)","","col");

  RunCanvas->cd(5);
  traw->Draw("adcList>>hADCList(196,0,196)","","");
  TPad *paddy5 = (TPad*)RunCanvas->FindObject("RunCanvas_5");
  paddy5->SetLogy();

  RunCanvas->cd(6);
  traw->Draw("tdcList>>hTDCList(128,0,128)","","");
  TPad *paddy6 = (TPad*)RunCanvas->FindObject("RunCanvas_6");
  paddy6->SetLogy();
}

void CheckAllRuns()
{
  OpenFile();

  TCanvas *Combo = new TCanvas("Combo","Combo",1800,1200);

  Combo->Divide(3,2);

  gROOT->ProcessLine(".x CUTposVdeltaE.C");
  TCutG *cuttySarc = (TCutG*)gROOT->FindObjectAny("CUTposVdeltaE");

  Combo->cd(1);
  trin->Draw("position:SiliconEnergy>>hCoinc(1000,0,4,512,0,4096)","SiliconEnergy<4 && TimeDifference>75 && TimeDifference<280","col");
  TPad *padMeHard1 = (TPad*)Combo->FindObject("Combo_1");
  padMeHard1->SetLogz();

  Combo->cd(2);
  trin->Draw("position:SiliconEnergy>>hCoincTritons(1000,0,4,512,0,4096)","SiliconEnergy<4 && TimeDifference>75 && TimeDifference<280 && CUTposVdeltaE","");
  TPad *padMeHard2 = (TPad*)Combo->FindObject("Combo_2");
  //padMeHard2->SetLogz();

  Combo->cd(3);
  trin->Draw("position:deltaE>>hPosdE(512,0,4096,512,0,4096)","","col");
  cuttySarc->Draw("same");

  Combo->cd(4);
  trin->Draw("16*(SiliconDetector-1)+SiliconPStrip:SiliconEnergy>>hSiliconEnergyPerStrip(1000,0,4,192,0,192)","TimeDifference>75 && TimeDifference<280","col");

  Combo->cd(5);
  int nentries = trin->GetEntries();
  char buffer[256];
  sprintf(buffer,"position:Entry$>>hPosVsEntryAll(%d,0,%d,512,0,4096)",nentries/100.,nentries);
  trin->Draw(buffer,"","col");

  Combo->cd(6);
  sprintf(buffer,"position:Entry$>>hPosVsEntryTritons(%d,0,%d,512,0,4096)",nentries/100.,nentries);
  //trin->Draw(buffer,"CUTposVdeltaE","col");

  gROOT->ProcessLine(".x CUTCoincLocust.C");
  trin->Draw("position>>hCoincEvents(512,0,4096)","TimeDifference>75 && TimeDifference<280 && CUTposVdeltaE && CUTCoincLocust","");

//  gROOT->ProcessLine(".x CutLocus1.C");
//  trin->Draw("position>>hCoincEvents(512,0,4096)","TimeDifference>75 && TimeDifference<280 && CUTposVdeltaE && CutLocus1","");
//  Combo->cd(1);
//  CutLocus1->Draw("same");
  

}

