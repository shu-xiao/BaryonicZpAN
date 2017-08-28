#include <iostream>
#include <algorithm>
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
    vector <vector<double>> punziList(2);
    vector <vector<double>> effi(2);
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
        effi[0].push_back(event[0][0]/nSigEvent);
        effi[1].push_back(event[1][0]/nSigEvent);
        punziList[0].push_back(punzi(event[0][0]/nSigEvent,event[0][1]));
        punziList[1].push_back(punzi(event[1][0]/nSigEvent,event[1][1]));
        //cout << effiF[i] << "\t" <<punziListF[i] << "\t" << event[0][1] << "\t" << TMath::Power(event[0][1],0.5)<< endl;
        //punziListF.push_back(SoverB(event[0][0],event[0][1]));
        //punziListR.push_back(SoverB(event[1][0],event[1][1]));
        //cout << "sig: " << event[1][0] << "\tbg: " << event[1][1] << "\tpz: " << punziListR[i] << endl;
    }
    reverse(effi[1].begin(),effi[1].end());
    reverse(punziList[1].begin(),punziList[1].end());
    for (int i=0;i<punziList[0].size();i++) {
        if (max[0]<punziList[0][i]) {
            max[0] = punziList[0][i];
            maxIndex[0] = i;
        }
        if (max[1]<punziList[1][i]) {
            max[1] = punziList[1][i];
            maxIndex[1] = i;
        }

    }
    return {binCentral, effi[0], effi[1], punziList[0], punziList[1]};
}
void setAxis(int color, float ymax, float x=1.0) {
    float xaxisPosition = (gPad->GetUxmax()-gPad->GetUxmin())*x+gPad->GetUxmin();
    TGaxis* axis = new TGaxis(xaxisPosition,gPad->GetUymin(),xaxisPosition,gPad->GetUymax(),0,ymax,510,"+L");
    axis->SetLineColor(color);
    axis->SetLabelColor(color);
    axis->Draw();
}
void setHist(TH1F* hist, int linecolor, int linewidth=2, bool addaxis=false, float xp=1.0) {
    Float_t rightmax = 1.1*hist->GetMaximum();
    Float_t scale    = gPad->GetUymax()/rightmax;
    hist->SetLineColor(linecolor);
    hist->SetLineWidth(linewidth);
    hist->Scale(scale);
    hist->Draw("histsame");
    if (addaxis) setAxis(linecolor,rightmax,xp);
    //hist->Scale(1/scale);
}


void ptOpt() {
    gStyle->SetStatY(0.9);
    gStyle->SetStatX(0.85);
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1","c1",800,600);
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
        h_bg[j] = (TH1F*)f_bg->Get(hname[j].Data());
        if (h_bg[j] == 0) cout << hname[j] << " cannot find matching figure in background root file" << endl;
        j++;
    }
    
    // draw for all figures
    TString pdfName = "optimize.pdf";
    vector <TH1F*> h_effi(2), h_punzi(2);
    TH1F *h_effiF, *h_effiR, *h_punziF, *h_punziR;
    vector<TString> hLogList = {"h_HT"};
    c1->Print((pdfName+"[").Data());
    for (int i=0; i<nHist;i++) {
        
        bool setlog = false;
        for (int ee=0;ee<hLogList.size();ee++) if (hname[i].Contains(hLogList[ee].Data())) setlog = true;
        c1->SetLogy(setlog);
        h_sig[i]->Scale(L2016/normL());
        h_sig[i]->SetLineColor(4);
        c1->Update();

        vector<vector<double>> info  = getCut(h_sig[i],h_bg[i]);
        int nBin = info[0].size();
        float binW = (info[0][1]-info[0][0]);
        h_effi[0] = new TH1F("h_effiLower","h_effiLower",nBin,info[0][0]-binW/2,info[0][nBin-1]+binW/2);
        h_effi[1] = new TH1F("h_effiHigher","h_effiHigher",nBin,info[0][0]-binW/2,info[0][nBin-1]+binW/2);
        h_punzi[0] = new TH1F("h_punziLower","h_punziLower",nBin,info[0][0]-binW/2,info[0][nBin-1]+binW/2);
        h_punzi[1] = new TH1F("h_punziHigher","h_punziHigher",nBin,info[0][0]-binW/2,info[0][nBin-1]+binW/2);
        for (int i=0;i<nBin;i++) {
            h_effi[0]->SetBinContent(i,info[1][i]);
            h_effi[1]->SetBinContent(i,info[2][i]);
            h_punzi[0]->SetBinContent(i,info[3][i]);
            h_punzi[1]->SetBinContent(i,info[4][i]);
        }
        
        cout << "i = " << i << "\thName = " << hname[i] << endl;
        vector<TH1F*> h_copy(2);
        TString suffix[2] = {"_lower","_higher"};
        for (int ee=0;ee<2;ee++) h_copy[ee] = (TH1F*)h_sig[i]->Clone("h_copy");
        for (j=0;j<2;j++) {
            c1->Clear();
            h_bg[i]->SetTitle((hname[i]+suffix[j]).Data());
            h_bg[i]->Draw("hist");
            c1->Update();
            setHist(h_punzi[j],51);
            setHist(h_effi[j],2,1,true);
            h_bg[i]->Draw("histsame");
            setHist(h_copy[j],4,3,true,0.9);
            
            if (i==0&&j==0) {
                TCanvas* c2 = new TCanvas("c2","c2",800,600);
                c2->cd();
                TLegend *legend = new TLegend(0.2,0.3,0.8,0.7);
                legend->SetHeader("The explanation for following figures","C");
                legend->AddEntry(h_sig[i],"signal");
                legend->AddEntry(h_bg[i],"QCD background");
                legend->AddEntry(h_effi[0],"efficiency vary with threshold");
                legend->AddEntry(h_punzi[0],"punzi value vary with threshold without axis");
                legend->Draw();
                c2->Print(pdfName.Data());
                delete c2;
            }
            c1->Print(pdfName.Data());
        }
    }
    c1->Print((pdfName+"]").Data());
}
