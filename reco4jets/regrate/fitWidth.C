#include <iostream>
#include "../setNCUStyle.C"

using namespace RooFit;
void fit(TH1F* h, int par=0) {
    string xTitle;
    static TLegend *le = new TLegend(0.15,0.7,0.45,0.9);
    float xMin, xMax;
    if (par==0) {
        xMin = -100;
        xMax = 100;
        le->SetX1NDC(0.2); 
        le->SetX2NDC(0.5); 
        xTitle = "M_{jj}-M_{h}";
    }
    if (par==1) {
        xMin = -300;
        xMax = 150;
        le->SetX1NDC(0.15); 
        le->SetX2NDC(0.45); 
        xTitle = "M_{jj}-M_{A0}";
    }
    if (par==2) {
        xMin = -400;
        xMax = 200;
        le->SetX1NDC(0.15); 
        le->SetX2NDC(0.45); 
        xTitle = "M_{jj}-M_{zp}";
    }
    le->Clear();
    RooRealVar x("x",xTitle.data(),0,xMin,xMax);
    RooRealVar x0_cb("x0_cb","x0_cb",0,xMin,xMax);
    RooRealVar x0_bw("x0_bw","x0_bw",0,xMin,xMax);
    RooDataHist hist_h("hist_h","hist_h",x,h);
    // CB
    RooRealVar n("n","n",-10,+10);
    RooRealVar sigma_cb("sigma","sigma",0,+200);
    RooRealVar sigma_bw("sigma_bw","sigma",0,+200);
    RooRealVar a("a","a",0,+200);
    RooCBShape cb("cb","cb",x,x0_cb,sigma_cb,a,n);
    cb.fitTo(hist_h);
    // BW
    RooBreitWigner bw("bw","bw",x,x0_bw,sigma_bw);
    bw.fitTo(hist_h);
    RooPlot* xframe = x.frame();
    hist_h.plotOn(xframe);
    cb.plotOn(xframe,LineColor(kBlue),RooFit::Name("CB"));
    bw.plotOn(xframe,LineColor(kRed),RooFit::Name("BW"));
    hist_h.plotOn(xframe);
    xframe->Draw();
    le->AddEntry(h,"genMatching","lp");
    le->AddEntry(xframe->findObject("CB"),"crystall ball","l");
    le->AddEntry(xframe->findObject("BW"),"BreitWigner","l");
    le->Draw();
}
void fitWidth() {
    setNCUStyle(1);
    gStyle->SetStatX(0.42);
    gStyle->SetStatY(0.9);
    gStyle->SetStatW(0.15);
    gStyle->SetStatFontSize(0.04);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(1);
    TCanvas *c1 = new TCanvas("c1","c1",4);
    TFile *f = TFile::Open("fitwidth.root");
    TH1F *h_hM = (TH1F*)f->Get("h_higgsM");
    TH1F *h_a0M = (TH1F*)f->Get("h_a0M");
    TH1F *h_zpM = (TH1F*)f->Get("h_ZpM");
    string name = "fitResult.pdf";
    c1->Print((name+"[").data());
    fit(h_hM,0);
    c1->Print(name.data());
    fit(h_a0M,1);
    c1->Print(name.data());
    fit(h_zpM,2);
    c1->Print(name.data());
    c1->Print((name+"]").data());
    //h_hM->Fit("gaus","","",-20,10);
    //h_hM->GetFunction("gaus")->SetLineWidth(3);
   
    
}
