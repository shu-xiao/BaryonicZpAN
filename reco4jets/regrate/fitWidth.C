#include <iostream>
#include "../setNCUStyle.C"

using namespace RooFit;
void fitWidth() {
    setNCUStyle(1);
    gStyle->SetStatX(0.42);
    gStyle->SetStatY(0.9);
    gStyle->SetStatW(0.15);
    gStyle->SetStatFontSize(0.04);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    TCanvas *c1 = new TCanvas("c1","c1",4);
    TFile *f = TFile::Open("fitwidth.root");
    TH1F *h_hM = (TH1F*)f->Get("h_higgsM");
    TH1F *h_a0M = (TH1F*)f->Get("h_a0M");
    TH1F *h_zpM = (TH1F*)f->Get("h_ZpM");
    string name = "fitResult.pdf";
    //c1->Print((name+"[").data());
    //h_hM->Fit("gaus","","",-20,10);
    //h_hM->GetFunction("gaus")->SetLineWidth(3);
    
    // using RooFit
    RooRealVar x("x","M_{jj}-M_{h}",0,-100,100);
    RooRealVar x0("x0","x0",0,-100,100);
    RooRealVar x0_bw("x0_bw","x0",0,-100,100);
    RooDataHist hist_h("hist_h","hist_h",x,h_hM);
    // CB
    RooRealVar n("n","n",-5,+5);
    RooRealVar sigma("sigma","sigma",0,+50);
    RooRealVar sigma_bw("sigma_bw","sigma",0,+50);
    RooRealVar a("a","a",0,+50);
    RooCBShape cb("cb","cb",x,x0,sigma,a,n);
    cb.fitTo(hist_h);
    // BW
    RooBreitWigner bw("bw","bw",x,x0_bw,sigma_bw);
    bw.fitTo(hist_h);
    // plot
    RooPlot* xframe = x.frame();
    xframe->SetName("");
    hist_h.plotOn(xframe);
    //hist_h.plotOn(xframe,DrawOption("hist"));
    cb.plotOn(xframe,LineColor(kBlue),RooFit::Name("CB"));
    bw.plotOn(xframe,LineColor(kRed),RooFit::Name("BW"));
    xframe->Draw();
    TLegend *le = new TLegend(0.2,0.7,0.5,0.9);
    le->AddEntry(h_hM,"genMatching","lp");
    le->AddEntry(xframe->findObject("CB"),"crystall ball","lp");
    le->AddEntry(xframe->findObject("BW"),"BreitWigner","lp");
    le->Draw();
    c1->Print(name.data());
    /*
    h_a0M->Fit("gaus","","",-70,10);
    h_a0M->GetFunction("gaus")->SetLineWidth(3);
    c1->Print(name.data());
    h_zpM->Fit("gaus","","",850,1000);
    h_zpM->GetFunction("gaus")->SetLineWidth(3);
    c1->Print(name.data());
    c1->Print((name+"]").data());
    */
    
}
