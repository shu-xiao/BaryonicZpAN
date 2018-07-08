#include <iostream>
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"

#define CISVV2_CUT 0.5426
using namespace std;
bool doReject = true;
void setMax(TH1F* h1, TH1F* h2) {
    float max = h1->GetMaximum();
    float max2 = h2->GetMaximum();
    if (max<max2) max = max2;
    h1->SetMaximum(max*1.05);
    h2->SetMaximum(max*1.05);
}
void setMax(TH1F* h1, TH1F* h2,TH1F* h3) {
    float max = h1->GetMaximum();
    float max2 = h2->GetMaximum();
    float max3 = h3->GetMaximum();
    if (max<max2) max = max2;
    if (max<max3) max = max3;
    h1->SetMaximum(max*1.4);
    h2->SetMaximum(max*1.4);
    h3->SetMaximum(max*1.4);
}
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
TH1F* copyHist(TH1F* h1, bool sig){
    TH1F* hclone = (TH1F*)h1->Clone(((string)h1->GetName()+Form("_copy%d",sig)).data());
    hclone->Reset();
    hclone->Sumw2();
    float bincenter, binError, binCon;
    for (int i=0;i<=h1->GetNbinsX();i++) {
        bincenter = h1->GetXaxis()->GetBinCenter(i);
        binError = h1->GetBinError(i);
        binCon = h1->GetBinContent(i);
        if (sig&&bincenter<140&&bincenter>120) {
            hclone->SetBinContent(i,binCon);
            hclone->SetBinError(i,binError);
        }
        else if (!sig&&(bincenter>140||bincenter<120)){
            hclone->SetBinContent(i,binCon);
            hclone->SetBinError(i,binError);

        }
    }
    if (sig) hclone->SetLineColor(kGreen);
    else hclone->SetLineColor(kBlack);
    hclone->SetLineWidth(2);
    return hclone;
}
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
void drawDiff(TH1F* h1, TH1F* h2, string title="",int io=0, string fileName="pfRatioCompare.pdf") {
    setMax(h1,h2);
    TCanvas* c3 = new TCanvas("c3","c3",3);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    c3->Divide(1,2,0.01,0.01);
    c3->cd(1);
    TH1F* h_copy = new TH1F(*h2);
    h_copy->Divide(h1);
    h2->SetLineColor(kBlack);
    h1->SetLineColor(kBlue);
    h1->GetYaxis()->SetTitle("A.U.");
    h1->GetXaxis()->SetTitle(title.data());
    c3->GetPad(1)->SetLeftMargin(0.12);
    c3->GetPad(2)->SetLeftMargin(0.12);
    c3->GetPad(1)->SetBottomMargin(0.12);
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->GetXaxis()->SetTitleSize(0.05);
    //h1->GetXaxis()->SetTitleOffset(0.4);
    //h1->GetYaxis()->SetTitleOffset(0.4);
    h2->SetStats(0);
    h1->Draw("e");
    h2->Draw("esame");
    // legend
    TLegend leg(0.5,0.7,0.88,0.85);
    TString htitle = h1->GetTitle();
    if (htitle.Contains("weight2")) leg.AddEntry(h1,"p/f(parabola funciton fit)*f");
    else leg.AddEntry(h1,"p/f(linear function fit)*f");
    leg.AddEntry(h2,"pass");
    leg.SetBorderSize(0);
    leg.Draw();
    // ratio
    c3->cd(2);
    c3->GetPad(2)->SetGridy();
    c3->GetPad(2)->SetPad(0.0,0.0,1,0.3);
    c3->GetPad(1)->SetPad(0.0,0.3,1,1);
    c3->GetPad(1)->SetTicks();
    c3->GetPad(2)->SetTicks();
    h_copy->GetYaxis()->SetRangeUser(0,2);
    h_copy->GetYaxis()->SetTitle("Ratio");
    h_copy->GetYaxis()->CenterTitle();
    h_copy->SetLineColor(kBlack);
    h_copy->GetXaxis()->SetLabelSize(0);
    h_copy->GetYaxis()->SetLabelSize(0.1);
    h_copy->GetYaxis()->SetTitleSize(0.1);
    h_copy->GetYaxis()->SetTitleOffset(0.5);
    h_copy->Draw("e");
    if (io==1) c3->Print((fileName+"[").data());
    c3->Print(fileName.data());
    if (io==2) c3->Print((fileName+"]").data());
    delete c3;
}
void treePF_btagandpt(bool doscan = false) {
    
    doscan = 1;
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1","c1",3);
    c1->SetLeftMargin(0.14);
    c1->SetRightMargin(0.14);
    TH2F* h2 = new TH2F("2f","QCD background",100,0,1,100,0,1000);
    h2->GetXaxis()->SetTitle("CISVV2");
    h2->GetYaxis()->SetTitle("Higgs Pt (GeV)");
    const float xsHTbeam[9] = {246400000,27990000,1712000,347700,32100,6831,1207,119.9,25.24};
    const float L2016=35.9*1000;//35.9 fb^-1
    const int maxCom = 50;
    string fName[] = {"tree_QCD_TMVA_HT50to100.root","tree_QCD_TMVA_HT100to200.root",\
    "tree_QCD_TMVA_HT200to300.root","tree_QCD_TMVA_HT300to500.root","tree_QCD_TMVA_HT500to700.root",\
    "tree_QCD_TMVA_HT700to1000.root","tree_QCD_TMVA_HT1000to1500.root",\
    "tree_QCD_TMVA_HT1500to2000.root","tree_QCD_TMVA_HT2000toInf.root"};
    string title[] = {"HT50to100","HT100to200","HT200to300","HT300to500","HT500to700",\
    "HT700to1000","HT1000to1500","HT1500to2000","HT2000toInf"};
    TH2F* hHT[9];
    for (int i=0;i<9;i++) hHT[i] = new TH2F(title[i].data(),title[i].data(),100,0,1,100,0,1000);
    Float_t CISVV2, CISVV2_1[maxCom], CISVV2_2[maxCom];
    Float_t hPt[maxCom],Mh[maxCom],hDeltaR[maxCom],hDeltaEta[maxCom], hDeltaPhi[maxCom],hptas[maxCom],hsdas[maxCom];
    float b = 1;
    int ij;
    string readFile;
    for (int ij=0;ij<9;ij++) {
        if (ij==0) continue;
        static Int_t nCom, nEvents;
        if (doscan) readFile = fName[ij];
        else readFile = "tree_0.root";
        TFile *f = new TFile(readFile.data());
        TH1F* h_allEvent = (TH1F*)f->Get("h_allEvent");
        nEvents = h_allEvent->GetEntries();
        if (doscan) b = L2016*xsHTbeam[ij]/nEvents;
        TTree *t1 = (TTree*)f->Get("tree");
        t1->SetBranchAddress("nCom",&nCom);
        t1->SetBranchAddress("CISVV2_Hb1",CISVV2_1);
        t1->SetBranchAddress("CISVV2_Hb2",CISVV2_2);
        t1->SetBranchAddress("Mh_weightLT",Mh);
        t1->SetBranchAddress("hPt",hPt);
        t1->SetBranchAddress("hDeltaR",hDeltaR);
        t1->SetBranchAddress("hDeltaPhi",hDeltaPhi);
        t1->SetBranchAddress("hDeltaEta",hDeltaEta);
        t1->SetBranchAddress("hptAs",hptas);
        t1->SetBranchAddress("hsdAs",hsdas);
        for (int i=0;i<t1->GetEntries();i++) {
            t1->GetEntry(i);
            bool isPass = false;
            int ind = 0;
            int find = -1;
            // search pass events
            for (int j=0;j<nCom;j++) {
                if (Mh[j]<90||Mh[j]>160) continue;
                if (CISVV2_1[j]>=CISVV2_CUT) {
                    isPass = true;
                    ind = j;
                    break;
                }
            }
            if (!isPass) {
                for (int j=0;j<nCom;j++) {
                    if (Mh[j]<90||Mh[j]>160) continue;
                    if (CISVV2_1[j]<CISVV2_CUT) {
                        find = j;
                        break;
                    }
                }
            }
            

            if (!isPass&&find<0) continue;
            h2->Fill(CISVV2_1[ind],hPt[ind],b);
            hHT[ij]->Fill(CISVV2_1[ind],hPt[ind],b);
        }
        f->Close();
        if (!doscan)break;
    } //end of fill
    c1->Print("hPttoBtag.pdf[");
    h2->Draw("colz");
    c1->Print("hPttoBtag.pdf");
    for (int i=0;i<9;i++) {
        hHT[i]->Draw("colz");
        c1->Print("hPttoBtag.pdf");
    }
    c1->Print("hPttoBtag.pdf]");

}

