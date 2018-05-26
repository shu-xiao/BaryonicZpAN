#include "TH1.h"
#include "TF1.h"
#include "TList.h"
#include <iostream>
#include "setNCUStyle.C"
using namespace std;
bool doReject = true;
double linear(double *x,double *par) {
    if (doReject&&x[0]<140&&x[0]>120) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0]+par[1]*x[0];
}
double pol_2(double *x,double *par) {
    if (doReject&&x[0]<140&&x[0]>120) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}
void pfRatioFit(string inputFile="QCD_TMVA_sample/QCDbg.root") {
    setNCUStyle(true);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    //gStyle->SetOptStat(0001112101.);
    gStyle->SetOptTitle(1);
    
    TCanvas *c1 = new TCanvas("c1","c1",3);
    TFile *f = TFile::Open(inputFile.data());
    TH1F *ratio = (TH1F*)f->Get("h_minhiggsM_ratio_L;2");
    TH1F *pass  = (TH1F*)f->Get("h_minhiggsM_pass_L;2");
    TH1F *fail  = (TH1F*)f->Get("h_minhiggsM_fail_L;2");
    TF1 *f1 = new TF1("f1",linear,40,200,2);
    TF1 *f2 = new TF1("f2",pol_2,40,200,3);
    f1->SetLineWidth(2);
    f2->SetLineWidth(2);
    f2->SetLineColor(kBlue);
    //f1->SetParameters(2,-1);
    ratio->Fit(f1);
    ratio->Fit(f2,"+");
    doReject = false;
    ratio->Draw();
    TLegend *leg = new TLegend(0.6,0.8,0.9,0.9);
    leg->AddEntry(f1,Form("linear, #chi^{2} = %.0f",f1->GetChisquare()),"l");
    leg->AddEntry(f2,Form("linear, #chi^{2} = %.0f",f2->GetChisquare()),"l");
    leg->Draw();
    /*
    doReject = false;
    // store 2 separate functions for visualization
    TF1 *f1eft = new TF1("f1eft",linier,0,120,2);
    f1eft->SetParameters(f1->GetParameters());
    ratio->GetListOfFunctions()->Add(f1eft);
    gROOT->GetListOfFunctions()->Remove(f1eft);
    TF1 *fright = new TF1("fright",linier,140,400,2);
    fright->SetParameters(f1->GetParameters());
    ratio->GetListOfFunctions()->Add(fright);
    gROOT->GetListOfFunctions()->Remove(fright);
    */
}
