#include <iostream>
#include "../setNCUStyle.C"
#include <RooBreitWigner.h>

//# define remote 0
#define largeRange
#define _DBCB 0
# ifdef remote
#undef _DBCB
#define _DBCB 1
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
# endif
using namespace RooFit;
void fit_(TH1F* h, int par=0) {
    string xTitle;
    static TLegend *le = new TLegend(0.15,0.75,0.45,0.9);
    float xMin, xMax;
#ifdef largeRange
    if (par==0) {
        xMin = -100;
        xMax = 100;
        le->SetX1NDC(0.2); 
        le->SetX2NDC(0.5); 
        xTitle = "M_{jj}-M_{h}";
    }
    if (par==1) {
        xMin = -300;
        xMax = 300;
        le->SetX1NDC(0.15); 
        le->SetX2NDC(0.45); 
        xTitle = "M_{jj}-M_{A0}";
    }
    if (par==2) {
        xMin = -500;
        xMax = 500;
        le->SetX1NDC(0.15); 
        le->SetX2NDC(0.45); 
        xTitle = "M_{jj}-M_{zp}";
    }
#else
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
#endif
    le->Clear();
    RooRealVar x("x",xTitle.data(),0,xMin,xMax);
    RooRealVar x0_cb("x0","x0",0,xMin,xMax);
    RooRealVar x0_bw("x0","x0",0,xMin,xMax);
    RooDataHist hist_h("hist_h","hist_h",x,h);
    // CB
    RooRealVar n("n","n",-10,+10);
    RooRealVar sigma_cb("sigma","sigma",0,+200);
    RooRealVar sigma_bw("sigma","sigma",0,+200);
    RooRealVar a("a","a",0,+200);
    RooCBShape cb("cb","cb",x,x0_cb,sigma_cb,a,n);
    cb.fitTo(hist_h);
    // BW
    RooBreitWigner bw("bw","bw",x,x0_bw,sigma_bw);
    bw.fitTo(hist_h);
    RooPlot* frame = x.frame(Name("frame"));
    hist_h.plotOn(frame);
    frame->SetMaximum(frame->GetMaximum()*1.1); 
    
    float fontSize = 0.025;
    cb.plotOn(frame,LineColor(kBlue),RooFit::Name("CB"));
    bw.plotOn(frame,LineColor(kRed),RooFit::Name("BW"));
    cb.paramOn(frame, Layout(0.65,0.93,0.73));
    frame->getAttText()->SetTextSize(fontSize);
    bw.paramOn(frame, Layout(0.65,0.93,0.48));
    frame->getAttText()->SetTextSize(fontSize);
    hist_h.plotOn(frame,Name("hist"));
    frame->Draw();
    // get chiSquare
    float x2_CB = frame->chiSquare("CB","hist",4);
    float x2_BW = frame->chiSquare("BW","hist",2);
    //float x2_CB = frame->chiSquare(4);
    //float x2_BW = frame->chiSquare(2);
    
    le->AddEntry(h,"genMatching","lp");
    le->AddEntry(frame->findObject("CB"),"crystall ball","l");
    le->AddEntry((TObject*)NULL,Form("#chi^{2} = %.3f",x2_CB),"");
    le->AddEntry(frame->findObject("BW"),"BreitWigner","l");
    le->AddEntry((TObject*)NULL,Form("#chi^{2} = %.3f",x2_BW),"");
    le->Draw();
}
# ifdef remote
void fitDBCBxBW(TH1F* h,int par=0) {
    string xTitle;
    static TLegend *le = new TLegend(0.15,0.75,0.40,0.9);
    le->Clear();
    float xMin, xMax;
#ifdef largeRange
    if (par==0) {
        xMin = -100;
        xMax = 100;
        le->SetX1NDC(0.2); 
        le->SetX2NDC(0.5); 
        xTitle = "M_{jj}-M_{h}";
    }
    if (par==1) {
        xMin = -300;
        xMax = 300;
        le->SetX1NDC(0.15); 
        le->SetX2NDC(0.45); 
        xTitle = "M_{jj}-M_{A0}";
    }
    if (par==2) {
        xMin = -500;
        xMax = 500;
        le->SetX1NDC(0.15); 
        le->SetX2NDC(0.4); 
        xTitle = "M_{jj}-M_{zp}";
    }
#else
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
#endif
    float mean = 0, sigma = 5;
    RooRealVar x("x",xTitle.data(),0, xMin, xMax);
    RooDataHist hist_h("hist_h","hist_h",x,h);
    //keep silence
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    RooMsgService::instance().setSilentMode(true);
      

    //BreitWigner
    RooRealVar Width_BW("#sigma_{BW}"," ",0,200);
    RooRealVar Mean_BW("#mu_{BW}"," ",0,xMin,xMax);
    RooBreitWigner BreitWigner("BreitWigner"," ",x,Mean_BW,Width_BW);

    //DBCB pdf var
    RooRealVar alpha1("#alpha_{1}"," ",0,200); // Alpha: Gaussian tail
    RooRealVar alpha2("#alpha_{2}"," ",0,200); // Alpha: Gaussian tail
    RooRealVar n1("n_{1}"," ",-10,+10);
    RooRealVar n2("n_{2}"," ",-10,+10);
    RooRealVar mean_CB("#mu_{DBCB}"," ",xMin,xMax);
    RooRealVar sigma_CB("#sigma_{DBCB}"," ",0 ,200);
    RooDoubleCB DBCB("x_pdf"," ", x, mean_CB, sigma_CB, alpha1, n1, alpha2, n2);
      
    //Convolution
    RooFFTConvPdf SigPdf("SigPdf"," ",x,BreitWigner,DBCB); // convolve
    SigPdf.fitTo(hist_h, Save(kTRUE), Range("x"),RooFit::NumCPU(1));
    
    //single DBCB pdf var
    RooRealVar alphaA("#alpha_{1}"," ",0,200); // Alpha: Gaussian tail
    RooRealVar alphaB("#alpha_{2}"," ",0,200); // Alpha: Gaussian tail
    RooRealVar nA("n_{1}"," ",-10,+10);
    RooRealVar nB("n_{2}"," ",-10,+10);
    RooRealVar mean_CB_A("#mu_{DBCB}"," ",xMin,xMax);
    RooRealVar sigma_CB_B("#sigma_{DBCB}"," ",0 ,200);
    RooDoubleCB DBCB_AB("x_pdf"," ", x, mean_CB_A, sigma_CB_B, alphaA, nA, alphaB, nB);
    DBCB_AB.fitTo(hist_h,Save(kTRUE),NumCPU(1));
    // fit and plot
    float fontSize = 0.022;
    RooPlot* frame = x.frame(Name("frame"));
    hist_h.plotOn(frame,Name("hist_h"));
    frame->SetMaximum(frame->GetMaximum()*1.1); 
    SigPdf.plotOn(frame,Name("curve"),LineColor(kBlue));
    DBCB_AB.plotOn(frame,Name("DBCB_AB"),LineColor(kRed));
    
    SigPdf.paramOn(frame,Layout(0.65,0.93,0.9));
    frame->getAttText()->SetTextSize(fontSize);
    DBCB_AB.paramOn(frame,Layout(0.15,0.36,0.74));
    frame->getAttText()->SetTextSize(fontSize);
    hist_h.plotOn(frame,Name("hist_h"));
    frame->Draw();

    float x2_DBCB = frame->chiSquare("curve","hist_h",8);
    float x2_DBCB_single = frame->chiSquare("DBCB_AB","hist_h",6);
    float fit_mean = mean_CB.getVal();
    float fit_sigma = sigma_CB.getVal();
    float fit_mean_Uncer = mean_CB.getError();
    float fit_sigma_Uncer = sigma_CB.getError();
    le->AddEntry(h,"genMatching","lp");
    le->AddEntry(frame->findObject("curve"),"DBCB#otimesBW","l");
    le->AddEntry((TObject*)0,Form("#chi^{2} = %.3f",x2_DBCB),"");
    le->AddEntry(frame->findObject("DBCB_AB"),"DBCB","l");
    le->AddEntry((TObject*)0,Form("#chi^{2} = %.3f",x2_DBCB_single),"");
    le->Draw();
}
#endif
void fitWidth() {
    setNCUStyle(1);
    gStyle->SetStatX(0.42);
    gStyle->SetStatY(0.9);
    gStyle->SetStatW(0.15);
    gStyle->SetStatFontSize(0.04);
    gStyle->SetOptFit(1111);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(1);
    TCanvas *c1 = new TCanvas("c1","c1",4);
    TFile *f = TFile::Open("fitwidth.root");
    TH1F *h_hM = (TH1F*)f->Get("h_higgsM");
    TH1F *h_a0M = (TH1F*)f->Get("h_a0M");
    TH1F *h_zpM = (TH1F*)f->Get("h_ZpM");
    
    void (*fit)(TH1F*,int);
    # ifdef remote
    fit = fitDBCBxBW; 
    string name = "fitResult_DBCB.pdf";
    # else 
    fit = fit_;
    string name = "fitResult.pdf";
    # endif
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
