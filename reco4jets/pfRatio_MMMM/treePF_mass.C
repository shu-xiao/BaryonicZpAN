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
    // h1 est, h2 sim
    setMax(h1,h2);
    TCanvas* c3 = new TCanvas("c3","c3",3);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    c3->Divide(1,2,0.01,0.01);
    c3->cd(1);
    TH1F* h_copy = new TH1F(*h1);
    h_copy->Divide(h2);
    h1->SetLineWidth(2);
    h2->SetLineWidth(2);
    h_copy->SetLineWidth(2);
    h1->SetLineColor(kBlue);
    h2->SetLineColor(kBlack);
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
    //if (htitle.Contains("weight2")) leg.AddEntry(h1,"p/f(parabola funciton fit)*f");
    //else leg.AddEntry(h1,"p/f(linear function fit)*f");
    // set estimate title
    if (htitle.Contains("weight2")) leg.AddEntry(h1,"estimation with parabola funciton");
    else leg.AddEntry(h1,"estimation with linear function");
    leg.AddEntry(h2,"simulation");
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
    h_copy->GetYaxis()->SetTitle("Estimation/MC");
    h_copy->GetYaxis()->CenterTitle();
    h_copy->SetLineColor(kBlack);
    h_copy->GetXaxis()->SetLabelSize(0);
    h_copy->GetYaxis()->SetLabelSize(0.1);
    h_copy->GetYaxis()->SetTitleSize(0.1);
    h_copy->GetYaxis()->SetTitleOffset(0.5);
    //h_copy->SetMarkerStyle(20);
    // Y = 1 LINE
    static TF1 *f1 = 0;
    if (!f1) f1 = new TF1("f1","1",-1000,1000);
    //f1->SetLineWidth(2);
    h_copy->Draw("E1");
    f1->Draw("same");
    //h_copy->Draw("e1same");
    if (io==1) c3->Print((fileName+"[").data());
    c3->Print(fileName.data());
    if (io==2) c3->Print((fileName+"]").data());
    delete c3;
}
//MainFun
void treePF_mass(bool doscan = false) {
    
    gStyle->SetOptStat(0);
    //doscan = 1;

    TCanvas *c1 = new TCanvas("c1","c1",3);
    TFile *ff = new TFile("ttt.root","recreate");
    const int nHist = 7;
    TH1F* h_Mh[nHist];
    TH1F* h_hPt[nHist];
    TH1F* h_hDeltaR[nHist];
    TH1F* h_hDeltaEta[nHist];
    TH1F* h_hDeltaPhi[nHist];
    TH1F* h_hptas[nHist], *h_hsdas[nHist];
    TH1F* h_MA0[nHist];
    TH1F* h_A0Pt[nHist];
    TH1F* h_A0DeltaR[nHist];
    TH1F* h_A0DeltaEta[nHist];
    TH1F* h_A0DeltaPhi[nHist];
    TH1F* h_A0ptas[nHist], *h_A0sdas[nHist];
    string suf[] = {"_fail","_pass","_testF","_testP","_ratio","_weight","_weight2"};
    for (int i=0;i<nHist;i++) {
        h_Mh[i] = new TH1F(Form("h_Mh%s",suf[i].data()),Form("h_Mh%s",suf[i].data()),20,90,160);
        h_hPt[i] = new TH1F(Form("h_hPt%s",suf[i].data()),Form("h_hPt%s",suf[i].data()),60,0,900);
        h_hDeltaR[i] = new TH1F(Form("h_hDeltaR%s",suf[i].data()),Form("h_hDeltaR%s",suf[i].data()),20,0,4);
        h_hDeltaPhi[i] = new TH1F(Form("h_hDeltaPhi%s",suf[i].data()),Form("h_hDeltaPhi%s",suf[i].data()),32,0,3.2);
        h_hDeltaEta[i] = new TH1F(Form("h_hDeltaEta%s",suf[i].data()),Form("h_hDeltaEta%s",suf[i].data()),20,0,4);
        h_hptas[i] = new TH1F(Form("h_hptas%s",suf[i].data()),Form("h_hptas%s",suf[i].data()),20,0,2);
        h_hsdas[i] = new TH1F(Form("h_hsdas%s",suf[i].data()),Form("h_hsdas%s",suf[i].data()),30,0,0.6);
        h_MA0[i] = new TH1F(Form("h_MA0%s",suf[i].data()),Form("h_MA0%s",suf[i].data()),20,250,350);
        h_A0Pt[i] = new TH1F(Form("h_A0Pt%s",suf[i].data()),Form("h_A0Pt%s",suf[i].data()),60,0,900);
        h_A0DeltaR[i] = new TH1F(Form("h_A0DeltaR%s",suf[i].data()),Form("h_A0DeltaR%s",suf[i].data()),20,0,4);
        h_A0DeltaPhi[i] = new TH1F(Form("h_A0DeltaPhi%s",suf[i].data()),Form("h_A0DeltaPhi%s",suf[i].data()),32,0,3.2);
        h_A0DeltaEta[i] = new TH1F(Form("h_A0DeltaEta%s",suf[i].data()),Form("h_A0DeltaEta%s",suf[i].data()),20,0,4);
        h_A0ptas[i] = new TH1F(Form("h_A0ptas%s",suf[i].data()),Form("h_A0ptas%s",suf[i].data()),20,0,2);
        h_A0sdas[i] = new TH1F(Form("h_A0sdas%s",suf[i].data()),Form("h_A0sdas%s",suf[i].data()),30,0,0.6);
        h_Mh[i]->Sumw2();
        h_hPt[i]->Sumw2();
        h_hDeltaR[i]->Sumw2();
        h_hDeltaPhi[i]->Sumw2();
        h_hDeltaEta[i]->Sumw2();
        h_hptas[i]->Sumw2();
        h_hsdas[i]->Sumw2();
        h_MA0[i]->Sumw2();
        h_A0Pt[i]->Sumw2();
        h_A0DeltaR[i]->Sumw2();
        h_A0DeltaPhi[i]->Sumw2();
        h_A0DeltaEta[i]->Sumw2();
        h_A0ptas[i]->Sumw2();
        h_A0sdas[i]->Sumw2();
    }
    TH1F* h_MA0_bin1 = new TH1F("h_MA0_bin1","h_MA0_bin1",20,250,350);
    TH1F* h_MA0_bin2 = new TH1F("h_MA0_bin2","h_MA0_bin2",20,250,350);
    
    const float xsHTbeam[9] = {246400000,27990000,1712000,347700,32100,6831,1207,119.9,25.24};
    const float L2016=35.9*1000;//35.9 fb^-1
    string fName[] = {"tree_QCD_TMVA_HT50to100.root","tree_QCD_TMVA_HT100to200.root","tree_QCD_TMVA_HT200to300.root","tree_QCD_TMVA_HT300to500.root","tree_QCD_TMVA_HT500to700.root","tree_QCD_TMVA_HT700to1000.root","tree_QCD_TMVA_HT1000to1500.root","tree_QCD_TMVA_HT1500to2000.root","tree_QCD_TMVA_HT2000toInf.root"};
    Float_t CISVV2_Hb1, CISVV2_Hb2, CISVV2_A0b1, CISVV2_A0b2;
    Float_t hPt,Mh,hDeltaR,hDeltaEta, hDeltaPhi,hptas,hsdas;
    Float_t A0Pt,MA0,A0DeltaR,A0DeltaEta, A0DeltaPhi,A0ptas,A0sdas;
    Bool_t isTag, isAntiTag;
    float b = 1;
    int ij;
    string readFile;
    for (int ij=0;ij<9;ij++) {
        if (ij==0&&1) continue;
        Int_t nCom, nEvents;
        if (doscan) readFile = fName[ij];
        else readFile = "../QCDsample_MMMM/QCD_HT100to200.root";
        TFile *f = new TFile(readFile.data());
        TH1F* h_allEvent = (TH1F*)f->Get("h_allEvent");
        nEvents = h_allEvent->GetEntries();
        if (doscan) b = L2016*xsHTbeam[ij]/nEvents;
        TTree *t1 = (TTree*)f->Get("tree");
        t1->SetBranchAddress("isTag",&isTag);
        t1->SetBranchAddress("isAntiTag",&isAntiTag);
        t1->SetBranchAddress("CISVV2_Hb1",&CISVV2_Hb1);
        t1->SetBranchAddress("CISVV2_Hb2",&CISVV2_Hb2);
        t1->SetBranchAddress("CISVV2_A0b1",&CISVV2_A0b1);
        t1->SetBranchAddress("CISVV2_A0b2",&CISVV2_A0b2);
        t1->SetBranchAddress("Mh",&Mh);
        t1->SetBranchAddress("hPt",&hPt);
        t1->SetBranchAddress("hDeltaR",&hDeltaR);
        t1->SetBranchAddress("hDeltaPhi",&hDeltaPhi);
        t1->SetBranchAddress("hDeltaEta",&hDeltaEta);
        t1->SetBranchAddress("hptAs",&hptas);
        t1->SetBranchAddress("hsdAs",&hsdas);
        t1->SetBranchAddress("MA0",&MA0);
        t1->SetBranchAddress("A0Pt",&A0Pt);
        t1->SetBranchAddress("A0DeltaR",&A0DeltaR);
        t1->SetBranchAddress("A0DeltaPhi",&A0DeltaPhi);
        t1->SetBranchAddress("A0DeltaEta",&A0DeltaEta);
        t1->SetBranchAddress("A0ptAs",&A0ptas);
        t1->SetBranchAddress("A0sdAs",&A0sdas);
        for (int i=0;i<t1->GetEntries();i++) {
            t1->GetEntry(i);
            bool isPass = isTag;
            // search pass events
            // if (Mh<90||Mh>160) continue;
            // search fail event
            h_Mh[i%2*2+isPass]->Fill(Mh,b);
            h_MA0[i%2*2+isPass]->Fill(MA0,b);
            h_hPt[i%2*2+isPass]->Fill(hPt,b);
            h_hDeltaR[i%2*2+isPass]->Fill(hDeltaR,b);
            h_hDeltaPhi[i%2*2+isPass]->Fill(hDeltaPhi,b);
            h_hDeltaEta[i%2*2+isPass]->Fill(hDeltaEta,b);
            h_hptas[i%2*2+isPass]->Fill(hptas,b);
            h_hsdas[i%2*2+isPass]->Fill(hsdas,b);
        }
        f->Close();
        if (!doscan)break;
    } //end of fill
    
    // ratio
    h_Mh[4]->Divide(h_Mh[1],h_Mh[0]);
    h_MA0[4]->Divide(h_MA0[1],h_MA0[0]);
    h_hPt[4]->Divide(h_hPt[1],h_hPt[0]);
    h_hDeltaR[4]->Divide(h_hDeltaR[1],h_hDeltaR[0]);
    h_hDeltaEta[4]->Divide(h_hDeltaEta[1],h_hDeltaEta[0]);
    h_hDeltaPhi[4]->Divide(h_hDeltaPhi[1],h_hDeltaPhi[0]);
    h_hptas[4]->Divide(h_hptas[1],h_hptas[0]);
    h_hsdas[4]->Divide(h_hsdas[1],h_hsdas[0]);
    ff->Write();
    // Fit

    TF1 *f1_s = new TF1("f1_s",linear,250,350,2);
    TF1 *f2_s = new TF1("f2_s",pol_2,250,350,3);
    f1_s->SetLineWidth(2);
    f2_s->SetLineWidth(2);
    f2_s->SetLineColor(kBlue);
    f1_s->SetLineColor(kRed);
    h_MA0[4]->Fit(f1_s);
    h_MA0[4]->Fit(f2_s,"+");
    string fileName = "treePFratio_mass.pdf";
    c1->Print((fileName+"[").data());
    doReject = false;
    h_MA0[4]->Draw("e");
    TLegend *leg = new TLegend(0.1,0.7,0.5,0.9);
    leg->AddEntry(h_MA0[4],"MA0 p/f ratio");
    leg->AddEntry(f1_s,Form("linear, #chi^{2}/ndf = %.1f/%d",f1_s->GetChisquare(),f1_s->GetNDF()),"l");
    leg->AddEntry(f2_s,Form("2nd order poly, #chi^{2}/ndf = %.1f/%d",f2_s->GetChisquare(),f2_s->GetNDF()),"l");
    leg->Draw();
    c1->Print(fileName.data());
    // empty SR
    leg->Clear();
    h_MA0[4]->Fit(f1_s);
    h_MA0[4]->Fit(f2_s,"+");
    TH1F* h_SR = copyHist(h_MA0[4],1);
    TH1F* h_SB = copyHist(h_MA0[4],0);
    h_SB->Draw("esame");
    h_SR->Draw("esame");
    leg->AddEntry(h_SB,"MA0 p/f ratio");
    leg->AddEntry(h_SR,"SR, not used in fitting");
    leg->AddEntry(f1_s,Form("linear, #chi^{2}/ndf = %.1f/%d",f1_s->GetChisquare(),f1_s->GetNDF()),"l");
    leg->AddEntry(f2_s,Form("2nd order poly, #chi^{2}/ndf = %.1f/%d",f2_s->GetChisquare(),f2_s->GetNDF()),"l");
    leg->Draw();
    c1->Print(fileName.data());
    h_SB->Draw("e");
    c1->Print(fileName.data());
    h_SR->Draw("e");
    c1->Print(fileName.data());

    // end
    h_MA0[0]->Draw("e");
    c1->Print(fileName.data());
    h_MA0[1]->Draw("e");
    c1->Print(fileName.data());
    // pf ratio
    h_MA0_bin1 = applyRatio(h_MA0[2],f1_s);
    h_MA0_bin2 = applyRatio(h_MA0[2],f2_s);
    h_MA0_bin1->SetLineColor(kRed);
    h_MA0[3]->SetLineColor(kBlue);
    h_MA0_bin2->SetLineColor(kCyan);
    setMax(h_MA0[3],h_MA0_bin1,h_MA0_bin2);
    h_MA0[3]->Draw("e");
    h_MA0_bin1->Draw("esame");
    h_MA0_bin2->Draw("esame");
    leg->Draw();
    c1->Print(fileName.data());
    
    // weight event by event
    b = 1;
    float bw, b2;
    for (int ij=0;ij<9;ij++) {
        if (ij==0&&1) continue;
        Int_t nCom, nEvents;
        if (doscan) readFile = fName[ij];
        else readFile = "../QCDsample_MMMM/QCD_HT100to200.root";
        TFile *f = new TFile(readFile.data());
        TH1F* h_allEvent = (TH1F*)f->Get("h_allEvent");
        nEvents = h_allEvent->GetEntries();
        if (doscan) b = L2016*xsHTbeam[ij]/nEvents;
        TTree *t1 = (TTree*)f->Get("tree");
        t1->SetBranchAddress("isTag",&isTag);
        t1->SetBranchAddress("isAntiTag",&isAntiTag);
        t1->SetBranchAddress("CISVV2_Hb1",&CISVV2_Hb1);
        t1->SetBranchAddress("CISVV2_Hb2",&CISVV2_Hb2);
        t1->SetBranchAddress("CISVV2_A0b1",&CISVV2_A0b1);
        t1->SetBranchAddress("CISVV2_A0b2",&CISVV2_A0b2);
        t1->SetBranchAddress("Mh",&Mh);
        t1->SetBranchAddress("hPt",&hPt);
        t1->SetBranchAddress("hDeltaR",&hDeltaR);
        t1->SetBranchAddress("hDeltaPhi",&hDeltaPhi);
        t1->SetBranchAddress("hDeltaEta",&hDeltaEta);
        t1->SetBranchAddress("hptAs",&hptas);
        t1->SetBranchAddress("hsdAs",&hsdas);
        t1->SetBranchAddress("MA0",&MA0);
        t1->SetBranchAddress("A0Pt",&A0Pt);
        t1->SetBranchAddress("A0DeltaR",&A0DeltaR);
        t1->SetBranchAddress("A0DeltaPhi",&A0DeltaPhi);
        t1->SetBranchAddress("A0DeltaEta",&A0DeltaEta);
        t1->SetBranchAddress("A0ptAs",&A0ptas);
        t1->SetBranchAddress("A0sdAs",&A0sdas);
        for (int i=0;i<t1->GetEntries();i++) {
            if (i%2==0) continue;
            t1->GetEntry(i);
            if (isTag) continue;
            if (MA0>350||MA0<250) continue;
            bool isPass = false;
            int find = -1;
            // veto pass, leave fail
            if (find<0) continue;
            bw = b*f1_s->Eval(MA0);
            b2 = b*f2_s->Eval(MA0);
            h_Mh[5]->Fill(Mh,bw);
            h_MA0[5]->Fill(MA0,bw);
            h_hPt[5]->Fill(hPt,bw);
            h_hDeltaR[5]->Fill(hDeltaR,bw);
            h_hDeltaPhi[5]->Fill(hDeltaPhi,bw);
            h_hDeltaEta[5]->Fill(hDeltaEta,bw);
            h_hptas[5]->Fill(hptas,bw);
            h_hsdas[5]->Fill(hsdas,bw);
            
            h_Mh[6]->Fill(Mh,b2);
            h_MA0[6]->Fill(MA0,b2);
            h_hPt[6]->Fill(hPt,b2);
            h_hDeltaR[6]->Fill(hDeltaR,b2);
            h_hDeltaPhi[6]->Fill(hDeltaPhi,b2);
            h_hDeltaEta[6]->Fill(hDeltaEta,b2);
            h_hptas[6]->Fill(hptas,b2);
            h_hsdas[6]->Fill(hsdas,b2);
        }
        f->Close();
        if (!doscan)break;
    } //end of fill
    h_Mh[5]->SetLineColor(kRed);
    h_MA0[5]->SetLineColor(kRed);
    h_hPt[5]->SetLineColor(kRed);
    h_hDeltaR[5]->SetLineColor(kRed);
    h_hDeltaEta[5]->SetLineColor(kRed);
    h_hDeltaPhi[5]->SetLineColor(kRed);
    h_hptas[5]->SetLineColor(kRed);
    h_hsdas[5]->SetLineColor(kRed);
    h_Mh[6]->SetLineColor(kCyan);
    h_MA0[6]->SetLineColor(kCyan);
    h_hPt[6]->SetLineColor(kCyan);
    h_hDeltaR[6]->SetLineColor(kCyan);
    h_hDeltaEta[6]->SetLineColor(kCyan);
    h_hDeltaPhi[6]->SetLineColor(kCyan);
    h_hptas[6]->SetLineColor(kCyan);
    h_hsdas[6]->SetLineColor(kCyan);
    
    leg->Clear();
    leg->SetX1NDC(0.45);
    leg->SetX2NDC(0.85);
    leg->SetY1NDC(0.75);
    leg->SetY2NDC(0.9);
    leg->AddEntry(h_hPt[3],"pass");
    leg->AddEntry(h_hPt[5],"fail*p/f (linier)");
    leg->AddEntry(h_hPt[6],"fail*p/f (2nd poly)");
    setMax(h_Mh[3],h_Mh[5],h_Mh[6]);
    h_Mh[3]->Draw("e");
    h_Mh[5]->Draw("esame");
    h_Mh[6]->Draw("esame");
    leg->Draw();
    c1->Print(fileName.data());
    setMax(h_MA0[3],h_MA0[5],h_MA0[6]);
    h_MA0[3]->Draw("e");
    h_MA0[5]->Draw("esame");
    h_MA0[6]->Draw("esame");
    leg->Draw();
    c1->Print(fileName.data());
    setMax(h_hPt[3],h_hPt[5],h_hPt[6]);
    h_hPt[3]->Draw("e");
    h_hPt[5]->Draw("esame");
    h_hPt[6]->Draw("esame");
    leg->Draw();
    c1->Print(fileName.data());
    setMax(h_hDeltaR[3],h_hDeltaR[5],h_hDeltaR[6]);
    h_hDeltaR[3]->Draw("e");
    h_hDeltaR[5]->Draw("esame");
    h_hDeltaR[6]->Draw("esame");
    leg->Draw();
    c1->Print(fileName.data());
    setMax(h_hDeltaEta[3],h_hDeltaEta[5],h_hDeltaEta[6]);
    h_hDeltaEta[3]->Draw("e");
    h_hDeltaEta[5]->Draw("esame");
    h_hDeltaEta[6]->Draw("esame");
    leg->Draw();
    c1->Print(fileName.data());
    setMax(h_hDeltaPhi[3],h_hDeltaPhi[5],h_hDeltaPhi[6]);
    h_hDeltaPhi[3]->Draw("e");
    h_hDeltaPhi[5]->Draw("esame");
    h_hDeltaPhi[6]->Draw("esame");
    leg->Draw();
    c1->Print(fileName.data());
    setMax(h_hptas[3],h_hptas[5],h_hptas[6]);
    h_hptas[3]->Draw("e");
    h_hptas[5]->Draw("esame");
    h_hptas[6]->Draw("esame");
    leg->Draw();
    c1->Print(fileName.data());
    setMax(h_hsdas[3],h_hsdas[5],h_hsdas[5]);
    h_hsdas[3]->Draw("e");
    h_hsdas[5]->Draw("esame");
    h_hsdas[6]->Draw("esame");
    leg->Draw();
    c1->Print(fileName.data());
    c1->Print((fileName+"]").data());
    // 3: sim,  5,6: estimate
    drawDiff(h_Mh[5],h_Mh[3],"M_{h}",1);
    drawDiff(h_Mh[6],h_Mh[3],"M_{h}");
    drawDiff(h_MA0[5],h_MA0[3],"M_{A0}");
    drawDiff(h_MA0[6],h_MA0[3],"M_{A0}");
    drawDiff(h_hPt[5],h_hPt[3],"higgs Pt");
    drawDiff(h_hPt[6],h_hPt[3],"higgs Pt");
    drawDiff(h_hDeltaR[5],h_hDeltaR[3],"#DeltaR_{bb}");
    drawDiff(h_hDeltaR[6],h_hDeltaR[3],"#DeltaR_{bb}");
    drawDiff(h_hDeltaEta[5],h_hDeltaEta[3], "#Delta#eta_{bb}");
    drawDiff(h_hDeltaEta[6],h_hDeltaEta[3],"#Delta#eta_{bb}");
    drawDiff(h_hDeltaPhi[5],h_hDeltaPhi[3],"#Delta#phi_{bb}");
    drawDiff(h_hDeltaPhi[6],h_hDeltaPhi[3],"#Delta#phi_{bb}");
    drawDiff(h_hptas[5],h_hptas[3],"pt assymetry (min(Pt1,Pt2)#Delta R/m_{jj})^{2}");
    drawDiff(h_hptas[6],h_hptas[3],"pt assymetry (min(Pt1,Pt2)#Delta R/m_{jj})^{2}");
    drawDiff(h_hsdas[5],h_hsdas[3],"min(pt1,pt2)/(pt1+pt2)");
    drawDiff(h_hsdas[6],h_hsdas[3],"min(pt1,pt2)/(pt1+pt2)",2);
}

