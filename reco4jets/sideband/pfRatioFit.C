#include "TH1.h"
#include "TF1.h"
#include "TList.h"
#include <iostream>
#include "setNCUStyle.C"
using namespace std;
bool doReject = true;
//vector <vector<float>> signalRange;
float signalRange[3][2]={{120,140},{},{}};
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
TH1F* divide(TH1F* hp, TH1F* hf, int nRebin = 1) {
    if (nRebin>1) {
        hp->Rebin(nRebin);
        hf->Rebin(nRebin);
    }
    TH1F* hclone = (TH1F*)hp->Clone(((string)hp->GetName()+"_ratio").data());
    hclone->Sumw2();
    hclone->Divide(hf);
    return hclone;
}
/*
TH1Fi* applyRatio(TH1F *hfail,TF1* fratio){
    TH1F* hclone = (TH1F*)hfail->Clone((hfail->GetName()+"_pfratio").data());
}
*/
TH1F* applyRatio(TH1F *hfail,TF1* fratio){
     TH1::SetDefaultSumw2(true);
    TH1F* hclone = (TH1F*)hfail->Clone(((string)hfail->GetName()+"_pfratio").data());
    hclone->Sumw2();
    float bincenter;
    for (int i=0;i<=hclone->GetNbinsX();i++){
        bincenter = hclone->GetXaxis()->GetBinCenter(i);
        //if (bincenter<140&&bincenter>120||true) {
            hclone->SetBinContent(i,hclone->GetBinContent(i)*fratio->Eval(bincenter));
        //}
    }
    return hclone;
}
void drawSig(TH1F* hist,float xmin=120,float xmax=140) {
    //hist->GetXaxis()->SetRangeUser(xmin,xmax);
    hist->SetAxisRange(120,140,"X");
    hist->SetLineColor(kGray);
    hist->Draw("esame");
}
void setrange() {
    //signalRange.push_back({120,140});
}
void pfRatioFit(string inputFile="QCD_TMVA_sample/QCDbg.root") {
    setNCUStyle(true);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    //gStyle->SetOptStat(0001112101.);
    //gStyle->SetOptTitle(1);
    static int k = 0;
    TCanvas *c1 = new TCanvas("c1","c1",3);
    TFile *f = TFile::Open(inputFile.data());
    TH1F *pass  = (TH1F*)f->Get("h_minhiggsM_pass_L");
    TH1F *fail  = (TH1F*)f->Get("h_minhiggsM_fail_L");
    TH1F *pass_copy = (TH1F*)pass->Clone("h_pass_copy");
    TH1F *fail_copy = (TH1F*)fail->Clone("h_fail_copy");
    pass_copy->SetAxisRange(100,160,"X");
    fail_copy->SetAxisRange(100,160,"X");
    //TH1F *ratio = (TH1F*)f->Get("h_minhiggsM_ratio_L;2");
    TH1F *ratio = divide(pass,fail,3);
    TH1F *ratio_s = divide(pass_copy,fail_copy);
    TF1 *f1 = new TF1("f1",linear,40,190,2);
    TF1 *f2 = new TF1("f2",pol_2,40,190,3);
    TF1 *f1_s = new TF1("f1_s",linear,100,160,2);
    TF1 *f2_s = new TF1("f2_s",pol_2,100,160,3);
    f1->SetLineWidth(2);
    f2->SetLineWidth(2);
    f1_s->SetLineWidth(2);
    f2_s->SetLineWidth(2);
    f2->SetLineColor(kBlue);
    f2_s->SetLineColor(kBlue);
    //f1->SetParameters(2,-1);
    ratio->Fit(f1);
    ratio->Fit(f2,"+");
    doReject = false;
    fail->Draw("e");
    c1->Print("ratioFit.pdf(");
    ratio->Draw();
    //drawSig(ratio);
    TLegend *leg = new TLegend(0.6,0.8,0.9,0.9);
    leg->AddEntry(ratio,"h_ratio");
    leg->AddEntry(f1,Form("linear, #chi^{2} = %.0f",f1->GetChisquare()),"l");
    leg->AddEntry(f2,Form("2nd order poly, #chi^{2} = %.0f",f2->GetChisquare()),"l");
    leg->Draw();
    c1->Print("ratioFit.pdf");
    // plot for short range
    doReject = true;
    ratio_s->GetXaxis()->SetTitle("M_{H}");
    ratio_s->Fit(f1_s);
    ratio_s->Fit(f2_s,"+");
    doReject = false;
    //fail_copy->Draw("e");
    //c1->Print("ratioFit.pdf");
    ratio_s->Draw();
    leg->Clear();
    leg->AddEntry(ratio_s,"h_ratio");
    leg->AddEntry(f1_s,Form("linear, #chi^{2} = %.2f",f1_s->GetChisquare()),"l");
    leg->AddEntry(f2_s,Form("2nd order poly, #chi^{2} = %.2f",f2_s->GetChisquare()),"l");
    leg->Draw();
    c1->Print("ratioFit.pdf");

    TH1F* h_ratio1 = applyRatio(fail,f1);
    TH1F* h_ratio2 = applyRatio(fail,f2);
    h_ratio1->SetLineColor(kBlue);
    h_ratio1->SetMarkerColor(kBlue);
    h_ratio2->SetLineColor(kRed);
    h_ratio2->SetMarkerColor(kRed);
    float ymax = max(h_ratio2->GetMaximum(),pass->GetMaximum());
    if (ymax>h_ratio1->GetMaximum()) h_ratio1->SetMaximum(ymax*1.3);
    h_ratio1->SetXTitle("M_{H}");
    h_ratio1->Draw("e");
    h_ratio2->Draw("esame");
    pass->Draw(" same");
    //drawSig(pass);
    leg->Clear();
    leg->AddEntry(pass,"h_pass");
    leg->AddEntry(h_ratio1,"linear");
    leg->AddEntry(h_ratio2,"poly_2");
    leg->Draw();
    c1->Print("ratioFit.pdf");
    
    TH1F* h_ratio1_s = applyRatio(fail_copy,f1_s);
    TH1F* h_ratio2_s = applyRatio(fail_copy,f2_s);
    h_ratio1_s->SetLineColor(kBlue);
    h_ratio1_s->SetMarkerColor(kBlue);
    h_ratio2_s->SetLineColor(kRed);
    h_ratio2_s->SetMarkerColor(kRed);
    ymax = max(h_ratio2_s->GetMaximum(),pass_copy->GetMaximum());
    if (ymax>h_ratio1_s->GetMaximum()) h_ratio1_s->SetMaximum(ymax*1.3);
    h_ratio1_s->SetXTitle("M_{H}");
    h_ratio1_s->Draw("e");
    h_ratio2_s->Draw("esame");
    pass_copy->Draw(" same");
    leg->Clear();
    leg->AddEntry(pass_copy,"h_pass");
    leg->AddEntry(h_ratio1_s,"linear");
    leg->AddEntry(h_ratio2_s,"poly_2");
    leg->Draw();
    c1->Print("ratioFit.pdf");
    c1->Clear();
    h_ratio1_s->DrawCopy("p");
    h_ratio2_s->DrawCopy("psame");
    pass_copy->DrawCopy("psame");
    c1->Update();
    leg->Draw();
    c1->Print("ratioFit.pdf)");
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
