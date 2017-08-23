#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TInterpreter.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#define L2016 35.9*1000 //35.9 fb^-1 = 35900 pb^-1
//#define L2016 1
using namespace std;
double punzi(double sigeff, double bg) { 
    return sigeff/(1+TMath::Power(bg,0.5));
}
double SoverB(double sig,double bg) {return sig/bg;}
float normL(double xs=0.0004199,int nEvent=10000) {
    // xs unit: pb
    return nEvent/xs; // pb^-1
}
vector <vector<double>> getCut(TH1F* h_sig, TH1F* h_bg) {
    double nSigEvent = h_sig->Integral();
    int nBin = h_sig->GetSize();
    vector <double> punziListF, punziListR;
    vector <double> effiF, effiR;
    vector <double> max(2);
    vector <double> binCentral;
    vector <int> maxIndex(2);
    double event[2][2] = {0};
    for (int i=0;i<nBin;i++) {
        binCentral.push_back(h_sig->GetBinCenter(i));
        event[0][0] += h_sig->GetBinContent(i);
        event[0][1] += h_bg->GetBinContent(i);
        event[1][0] += h_sig->GetBinContent(nBin-i-1);
        event[1][1] += h_bg->GetBinContent(nBin-i-1);
        effiF.push_back(event[0][0]/nSigEvent);
        effiR.push_back(event[1][0]/nSigEvent);
        punziListF.push_back(punzi(event[0][0]/nSigEvent,event[0][1]));
        punziListR.push_back(punzi(event[1][0]/nSigEvent,event[1][1]));
        //cout << effiF[i] << "\t" <<punziListF[i] << "\t" << event[0][1] << "\t" << TMath::Power(event[0][1],0.5)<< endl;
        //punziListF.push_back(SoverB(event[0][0],event[0][1]));
        //punziListR.push_back(SoverB(event[1][0],event[1][1]));
        //cout << "sig: " << event[1][0] << "\tbg: " << event[1][1] << "\tpz: " << punziListR[i] << endl;
    }
    for (int i=0;i<punziListF.size();i++) {
        if (max[0]<punziListF[i]) {
            max[0] = punziListF[i];
            maxIndex[0] = i;
        }
        if (max[1]<punziListR[i]) {
            max[1] = punziListR[i];
            maxIndex[1] = i;
        }

    }
    return {binCentral, effiF, effiR, punziListF, punziListR};
}
void gpd() {cout << gPad->GetUymax() << "te"<< endl; }
void ptOpt() {
    gStyle->SetStatY(0.9);
    gStyle->SetStatX(0.85);
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1","c1",800,600);
    //TCanvas* c2 = new TCanvas("c2","c2",800,600);
    //c1->cd();
    TFile* f_signal = TFile::Open("../signal/bb2HDM_MZp800_MA0300.root");
    TFile* f_bg = TFile::Open("../QCDsample/QCDbg.root");
    
    // get histogram name list

    const int nHist = 25; 
    TH1F *h_sig[nHist], *h_bg[nHist];
    TString hname[nHist];
    TIter keyList(f_signal->GetListOfKeys()); 
    TKey *key;
    int j = 0;
    while ((key = (TKey*)keyList())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1")) continue;
        h_sig[j] = (TH1F*)key->ReadObj();
        if (h_sig[j]==0) cout << hname[j] << "cannot be loaded from signal root file" << endl;
        hname[j] = h_sig[j]->GetName();
        if (hname[j].Contains("h_ZptoHA0DeltaR")) hname[j] = "h_zptobbDeltaR";
        if (hname[j].Contains("h_HiggstobbDeltaEta")) hname[j] = "h_HiggstobbdeltaEta";
        if (hname[j].Contains("h_HiggstobbDeltaPhi")) hname[j] = "h_HiggstobbdeltaPhi";
        if (hname[j].Contains("h_A0tobbDeltaPhi")) hname[j] = "h_A0tobbdeltaPhi";
        if (hname[j].Contains("h_A0tobbDeltaEta")) hname[j] = "h_A0tobbdeltaEta";
        if (hname[j].Contains("h_ZptoHA0DeltaEta")) hname[j] = "h_ZptobbdeltaEta";
        if (hname[j].Contains("h_ZptoHA0DeltaPhi")) hname[j] = "h_ZptobbdeltaPhi";
        h_bg[j] = (TH1F*)f_bg->Get(hname[j].Data());
        if (h_bg[j] == 0) cout << hname[j] << " cannot find matching figure in background root file" << endl;
        j++;
    }
    
    // draw for all figures
    TString pdfName = "optimize.pdf";
    c1->Print((pdfName+"[").Data());
    TH1F *h_effiF, *h_effiR, *h_punziF, *h_punziR;
    vector<TString> hLogList = {"h_HT"};
    for (int i=0; i<nHist;i++) {
        
        c1->Clear();
        bool setlog = false;
        for (int ee=0;ee<hLogList.size();ee++) if (hname[i].Contains(hLogList[ee].Data())) setlog = true;
        c1->SetLogy(setlog);
        h_sig[i]->Scale(L2016/normL());
        h_sig[i]->SetLineColor(4);
        h_bg[i]->Draw("hist");
        c1->Update();

        vector<vector<double>> info  = getCut(h_sig[i],h_bg[i]);
        int nBin = info[0].size();
        float binW = (info[0][1]-info[0][0]);
        h_effiF = new TH1F("h_effiF","h_effiF",nBin,info[0][0]-binW/2,info[0][nBin-1]+binW/2);
        h_effiR = new TH1F("h_effiR","h_effiR",nBin,info[0][0]-binW/2,info[0][nBin-1]+binW/2);
        h_punziF = new TH1F("h_punziF","h_punziF",nBin,info[0][0]-binW/2,info[0][nBin-1]+binW/2);
        h_punziR = new TH1F("h_punziR","h_punziR",nBin,info[0][0]-binW/2,info[0][nBin-1]+binW/2);
        for (int i=0;i<nBin;i++) {
            h_effiF->SetBinContent(i,info[1][i]);
            h_effiR->SetBinContent(i,info[2][i]);
            h_punziF->SetBinContent(i,info[3][i]);
            h_punziR->SetBinContent(i,info[4][i]);
        }
        c1->Update();
        
        Float_t rightmax3 = 1.1*h_punziF->GetMaximum();
        Float_t scale3    = gPad->GetUymax()/rightmax3;
        h_punziF->SetLineColor(51);
        h_punziF->SetLineWidth(2);
        h_punziF->Scale(scale3);
        h_punziF->Draw("histsame");
        c1->Update();
        
        Float_t rightmax = 1.1*h_effiF->GetMaximum();
        Float_t scale    = gPad->GetUymax()/rightmax;
        h_effiF->SetLineColor(kRed);
        h_effiF->SetLineWidth(1);
        h_effiF->Scale(scale);
        h_effiF->Draw("histsame");
        //draw an axis on the right side
        TGaxis* axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),0,rightmax,510,"+L");
        axis->SetLineColor(kRed);
        axis->SetLabelColor(kRed);
        
        
        h_bg[i]->Draw("histsame");

        int lineColor2 = 4;
        Float_t rightmax2 = 1.1*h_sig[i]->GetMaximum();
        Float_t scale2    = gPad->GetUymax()/rightmax2;
        h_sig[i]->SetLineColor(lineColor2);
        h_sig[i]->SetLineWidth(3);
        h_sig[i]->Scale(scale2);
        h_sig[i]->Draw("histsame");
        
        //draw an axis on the right side
        float xaxisPosition = (gPad->GetUxmax()-gPad->GetUxmin())*0.9+gPad->GetUxmin();
        TGaxis* axis2 = new TGaxis(xaxisPosition,gPad->GetUymin(),xaxisPosition,gPad->GetUymax(),0,rightmax2,510,"+L");
        axis2->SetLineColor(lineColor2);
        axis2->SetLabelColor(lineColor2);
        axis->Draw();
        axis2->Draw();
        cout << "i = " << i << "\thName = " << hname[i] << endl;
        
        if (i==0) {
            TCanvas* c2 = new TCanvas("c2","c2",800,600);
            c2->cd();
            TLegend *legend = new TLegend(0.2,0.3,0.8,0.7);
            legend->SetHeader("The explanation for following figures","C");
            legend->AddEntry(h_sig[i],"signal");
            legend->AddEntry(h_bg[i],"QCD background");
            legend->AddEntry(h_effiF,"efficiency vary with threshold");
            legend->AddEntry(h_punziF,"punzi value vary with threshold without axis");
            legend->Draw();
            c2->Print(pdfName.Data());
            delete c2;
        }
        gpd();
        c1->Print(pdfName.Data());
    }
    c1->Print((pdfName+"]").Data());
}
