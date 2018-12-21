#include <iostream>
//#include "RooParametricHist2D.h"
using namespace RooFit;
void readWorkspace() {
    TFile *f = new TFile("base_4bMMMM_nosig.root");
    RooWorkspace* w = (RooWorkspace*)f->Get("w_2D");
    RooRealVar* ptrMh = w->var("Mh") ;
    RooRealVar* ptrMA0 = w->var("MA0") ;
    RooRealVar  MA0;
    RooRealVar  Mh("Mh","Mh",0);
    //MA0 = *ptrMA0;
    //Mh = *ptrMh;
    
    RooAbsPdf* qcd_fail = w->pdf("qcd_fail") ;
    RooAbsPdf* qcd_pass = w->pdf("qcd_pass") ;
    cout << "test" << endl;
    RooDataHist* data_obs_fail = (RooDataHist*)w->data("data_obs_fail") ;
    RooDataHist* data_obs_pass = (RooDataHist*)w->data("data_obs_pass") ;
    //RooDataSet* data_obs_fail = (RooDataSet*)w->data("data_obs_fail") ;
    //RooDataSet* data_obs_pass = (RooDataSet*)w->data("data_obs_pass") ;

    //TH2F* h_2dpass = data_obs_pass->createHistogram(MA0,Mh,"","h_2dpass");
    //TH2F* h_2dfail = data_obs_fail->createHistogram(MA0,Mh,"","h_2dfail");
    //TH2F* h_2dpass = data_obs_pass->createHistogram(*MA0,*Mh,"","h_2dpass");
    //TH2F* h_2dfail = data_obs_fail->createHistogram(*MA0,*Mh,"","h_2dfail");
    //TH1* h_2dpass = data_obs_pass->createHistogram("h_2dpass",MA0,Binning(20),YVar(Mh,Binning(30)));
}
