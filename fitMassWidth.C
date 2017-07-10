#include "TH1.h"
#include "TF1.h"
#include "TList.h"
#include "TMathBase.h"
#include "TMath.h"
#include <TFile.h>
#include "TROOT.h"
#include <TCanvas.h>
#include <string>
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"

#include "RooPlot.h"

#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "RooTFnBinding.h" 
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "TVirtualFFT.h"
#include "RooBifurGauss.h"
#include "RooHistPdf.h"
#include "RooBreitWigner.h"
#include "RooNovosibirsk.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "TText.h"
#include "RooChi2Var.h"
#include "RooCBShape.h"
#include "RooVoigtian.h"

using namespace RooFit;
void fitMassWidth(std::string rootFile) {
    TCanvas *c1 = new TCanvas("c1","",1440,900);
    TFile *f1 = TFile::Open(rootFile.data());
    TH1F *h_zpMass = (TH1F*) f1->Get("h_ZpJetM");
    RooRealVar x("x","x",0,2500);
    x.setRange("r1",0,1400);
    x.setRange("r2",1500,2500);
    x.setRange("r3",1200,180);
    RooDataHist dh("dh","dh",x,Import(*h_zpMass)) ;
    RooPlot* frame;
    frame = x.frame(Title("gauss"));
    dh.plotOn(frame);
    
    RooRealVar ml("ml","mean landau",150,1,200) ;
    RooRealVar sl("sl","sigma landau",40,0.1,100) ;
  
    RooRealVar mean("mean","mean",240,0,2500) ;
    RooRealVar sigma("sigma","sigma",200,0.1,500) ;
    RooGaussian gauss("gauss","gauss",x,mean,sigma) ;

    RooRealVar meanW("meanW", "meanW",100,0,2500) ;
    RooRealVar width("width","width",20,0.1,400) ;
    RooBreitWigner BW("BW","BW",x,meanW,width);
    
    RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.1,0.,1.) ;
    RooLandau landauA("landauA","landauA",x,ml,sl) ;
    RooFFTConvPdf lxgB("lxgB","landau (X) BreitWigner",x,gauss,BW);
    RooAddPdf Exp2("Exp2","Exp2",RooArgList(lxgB,landauA),sig1frac); 
    lxgB.fitTo(dh,Range("r3"));
    landauA.fitTo(dh,Range("r1","r2"));

    Exp2.fitTo(dh);
    Exp2.plotOn(frame);
    Exp2.paramOn(frame,Layout(0.55));
    Exp2.plotOn(frame,Components("lxgB")) ;
    Exp2.plotOn(frame,Components(landauA),LineStyle(kDashed)) ;
    frame->Draw();
}
