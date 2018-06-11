#include <iostream>
#define CISVV2_CUT 0.5426
using namespace std;
bool doReject = true;
void setMax(TH1F* h1, TH1F* h2) {
    float max = h1->GetMaximum();
    float max2 = h2->GetMaximum();
    if (max<max2) max = max2;
    h1->SetMaximum(max*1.1);
    h2->SetMaximum(max*1.1);
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
void treePF_mass(bool doscan = false) {
    
    TCanvas *c1 = new TCanvas("c1","c1",3);
    TFile *ff = new TFile("ttt.root","recreate");
    const int nHist = 7;
    TH1F* h_Mh[nHist];
    TH1F* h_hPt[nHist];
    TH1F* h_hDeltaR[nHist];
    TH1F* h_hDeltaEta[nHist];
    TH1F* h_hDeltaPhi[nHist];
    TH1F* h_hptas[nHist], *h_hsdas[nHist];
    string suf[] = {"_fail","_pass","_testF","_testP","_ratio","_weight","_weight2"};
    for (int i=0;i<nHist;i++) {
        h_Mh[i] = new TH1F(Form("h_Mh%s",suf[i].data()),Form("h_Mh%s",suf[i].data()),20,90,160);
        h_hPt[i] = new TH1F(Form("h_hPt%s",suf[i].data()),Form("h_hPt%s",suf[i].data()),60,0,900);
        h_hDeltaR[i] = new TH1F(Form("h_hDeltaR%s",suf[i].data()),Form("h_hDeltaR%s",suf[i].data()),20,0,4);
        h_hDeltaPhi[i] = new TH1F(Form("h_hDeltaPhi%s",suf[i].data()),Form("h_hDeltaPhi%s",suf[i].data()),32,0,3.2);
        h_hDeltaEta[i] = new TH1F(Form("h_hDeltaEta%s",suf[i].data()),Form("h_hDeltaEta%s",suf[i].data()),20,0,4);
        h_hptas[i] = new TH1F(Form("h_hptas%s",suf[i].data()),Form("h_hptas%s",suf[i].data()),20,0,2);
        h_hsdas[i] = new TH1F(Form("h_hsdas%s",suf[i].data()),Form("h_hsdas%s",suf[i].data()),30,0,0.6);
        h_Mh[i]->Sumw2();
        h_hPt[i]->Sumw2();
        h_hDeltaR[i]->Sumw2();
        h_hDeltaPhi[i]->Sumw2();
        h_hDeltaEta[i]->Sumw2();
        h_hptas[i]->Sumw2();
        h_hsdas[i]->Sumw2();
    }

    const float xsHTbeam[9] = {246400000,27990000,1712000,347700,32100,6831,1207,119.9,25.24};
    const float L2016=35.9*1000;//35.9 fb^-1
    const int maxCom = 50;
    string fName[] = {"tree_QCD_TMVA_HT50to100.root","tree_QCD_TMVA_HT100to200.root","tree_QCD_TMVA_HT200to300.root","tree_QCD_TMVA_HT300to500.root","tree_QCD_TMVA_HT500to700.root","tree_QCD_TMVA_HT700to1000.root","tree_QCD_TMVA_HT1000to1500.root","tree_QCD_TMVA_HT1500to2000.root","tree_QCD_TMVA_HT2000toInf.root"};
    Float_t CISVV2, CISVV2_1[maxCom], CISVV2_2[maxCom];
    Float_t hPt[maxCom],Mh[maxCom],hDeltaR[maxCom],hDeltaEta[maxCom], hDeltaPhi[maxCom],hptas[maxCom],hsdas[maxCom];
    float b = 1;
    int ij;
    string readFile;
    for (int ij=0;ij<9;ij++) {
        if (ij==0&&1) continue;
        Int_t nCom, nEvents;
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
            if (isPass) {
                h_Mh[i%2*2+isPass]->Fill(Mh[ind],b);
                h_hPt[i%2*2+isPass]->Fill(hPt[ind],b);
                h_hDeltaR[i%2*2+isPass]->Fill(hDeltaR[ind],b);
                h_hDeltaPhi[i%2*2+isPass]->Fill(hDeltaPhi[ind],b);
                h_hDeltaEta[i%2*2+isPass]->Fill(hDeltaEta[ind],b);
                h_hptas[i%2*2+isPass]->Fill(hptas[ind],b);
                h_hsdas[i%2*2+isPass]->Fill(hsdas[ind],b);
            }
            else {
                if (find<0) continue;
                h_Mh[i%2*2+isPass]->Fill(Mh[find],b);
                h_hPt[i%2*2+isPass]->Fill(hPt[find],b);
                h_hDeltaR[i%2*2+isPass]->Fill(hDeltaR[find],b);
                h_hDeltaPhi[i%2*2+isPass]->Fill(hDeltaPhi[find],b);
                h_hDeltaEta[i%2*2+isPass]->Fill(hDeltaEta[find],b);
                h_hptas[i%2*2+isPass]->Fill(hptas[find],b);
                h_hsdas[i%2*2+isPass]->Fill(hsdas[find],b);
            }
        }
        f->Close();
        if (!doscan)break;
    } //end of fill
    
    // ratio
    h_Mh[4]->Divide(h_Mh[1],h_Mh[0]);
    h_hPt[4]->Divide(h_hPt[1],h_hPt[0]);
    h_hDeltaR[4]->Divide(h_hDeltaR[1],h_hDeltaR[0]);
    h_hDeltaEta[4]->Divide(h_hDeltaEta[1],h_hDeltaEta[0]);
    h_hDeltaPhi[4]->Divide(h_hDeltaPhi[1],h_hDeltaPhi[0]);
    h_hptas[4]->Divide(h_hptas[1],h_hptas[0]);
    h_hsdas[4]->Divide(h_hsdas[1],h_hsdas[0]);
    ff->Write();
    // Fit

    TF1 *f1_s = new TF1("f1_s",linear,100,160,2);
    TF1 *f2_s = new TF1("f2_s",pol_2,100,160,3);
    f1_s->SetLineWidth(2);
    f2_s->SetLineWidth(2);
    f2_s->SetLineColor(kBlue);
    f1_s->SetLineColor(kRed);
    h_Mh[4]->Fit(f1_s);
    h_Mh[4]->Fit(f2_s,"+");
    string fileName = "treePFratio_mass.pdf";
    c1->Print((fileName+"[").data());
    doReject = false;
    h_Mh[4]->Draw("e");
    TLegend *leg = new TLegend(0.1,0.7,0.5,0.9);
    leg->AddEntry(h_Mh[4],"Mh p/f ratio");
    leg->AddEntry(f1_s,Form("linear, #chi^{2}/ndf = %.1f/%d",f1_s->GetChisquare(),f1_s->GetNDF()),"l");
    leg->AddEntry(f2_s,Form("2nd order poly, #chi^{2}/ndf = %.1f/%d",f2_s->GetChisquare(),f2_s->GetNDF()),"l");
    leg->Draw();
    c1->Print(fileName.data());
    h_Mh[5] = applyRatio(h_Mh[2],f1_s);
    h_Mh[6] = applyRatio(h_Mh[2],f2_s);
    h_Mh[5]->SetLineColor(kRed);
    h_Mh[3]->SetLineColor(kBlue);
    h_Mh[6]->SetLineColor(kCyan);
    setMax(h_Mh[3],h_Mh[5],h_Mh[6]);
    h_Mh[3]->Draw("e");
    h_Mh[5]->Draw("esame");
    h_Mh[6]->Draw("esame");
    leg->Draw();
    c1->Print(fileName.data());
    
    // weight event by event
    b = 1;
    float bw, b2;
    for (int ij=0;ij<9;ij++) {
        if (ij==0&&1) continue;
        Int_t nCom, nEvents;
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
            if (!i%2) continue;
            t1->GetEntry(i);
            bool isPass = false;
            int find = -1;
            for (int j=0;j<nCom;j++) {
                
                if (Mh[j]<90||Mh[j]>160) continue;
                if (CISVV2_1[j]>=CISVV2_CUT) {
                    isPass = true;
                    break;
                }
            }
            if (isPass) continue;
            for (int j=0;j<nCom;j++) {
                if (Mh[j]<90||Mh[j]>160) continue;
                if (CISVV2_1[j]<CISVV2_CUT) {
                    find = j;
                    break;
                }
            }
            if (find<0) continue;
            bw = b*f1_s->Eval(Mh[find]);
            b2 = b*f2_s->Eval(Mh[find]);
            h_hPt[5]->Fill(hPt[find],bw);
            h_hDeltaR[5]->Fill(hDeltaR[find],bw);
            h_hDeltaPhi[5]->Fill(hDeltaPhi[find],bw);
            h_hDeltaEta[5]->Fill(hDeltaEta[find],bw);
            h_hptas[5]->Fill(hptas[find],bw);
            h_hsdas[5]->Fill(hsdas[find],bw);
            
            h_hPt[6]->Fill(hPt[find],b2);
            h_hDeltaR[6]->Fill(hDeltaR[find],b2);
            h_hDeltaPhi[6]->Fill(hDeltaPhi[find],b2);
            h_hDeltaEta[6]->Fill(hDeltaEta[find],b2);
            h_hptas[6]->Fill(hptas[find],b2);
            h_hsdas[6]->Fill(hsdas[find],b2);
        }
        f->Close();
        if (!doscan)break;
    } //end of fill
    h_hPt[5]->SetLineColor(kRed);
    h_hDeltaR[5]->SetLineColor(kRed);
    h_hDeltaEta[5]->SetLineColor(kRed);
    h_hDeltaPhi[5]->SetLineColor(kRed);
    h_hptas[5]->SetLineColor(kRed);
    h_hsdas[5]->SetLineColor(kRed);
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
}
