#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TInterpreter.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
//#define L2016 35.9*1000 //35.9 pb^-1
#define L2016 1
using namespace std;
double punzi(double sigeff, double bg) { 
    return sigeff/(1+TMath::Power(bg,1/2));
}
double SoverB(double sig,double bg) {return sig/bg;}
float normL(double xs=0.0004199,int nEvent=10000) {
    // xs unit: pb
    return nEvent/xs;
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
    /*
    TGraph *tg = new TGraph(nBin,&binCentral[0],&punziListF[0]);
    TCanvas* c2 = new TCanvas("c2","c2",800,600);
    tg->Draw("AP*");
    c2->SaveAs("pz.png");
    */
    //return {h_sig->GetBinCenter(maxIndex[0]),h_sig->GetBinCenter(maxIndex[1])};
    return {binCentral, effiF, effiR, punziListF, punziListR};
}
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

    const int nHist = 19; 
    TH1F *h_sig[nHist], *h_bg[nHist];
    TString hname[nHist];
    TIter keyList(f_signal->GetListOfKeys()); 
    TKey *key;
    int j = 0;
    while ((key = (TKey*)keyList())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1")) continue;
        h_sig[j] = (TH1F*)key->ReadObj();
        hname[j] = h_sig[j]->GetName();
        if (hname[j].Contains("h_HT")) {
            continue;
        }
        h_bg[j] = (TH1F*)f_bg->Get(hname[j].Data());
        j++;
    }
    
    // draw for all figures
    TString pdfName = "optimize.pdf";
    c1->Print((pdfName+"[").Data());
    TH1F *h_effiF, *h_effiR, *h_punziF, *h_punziR;
    for (int i=0; i<nHist;i++) {
        
        
        c1->Clear();
        //TH1F* h_ptsig = (TH1F*) f_signal->Get("h_higgsPt");    
        //TH1F* h_ptbg = (TH1F*) f_bg->Get("h_higgsPt");
        //h_ptsig->Scale(L2016/normL());
        //h_ptsig->Scale(1/h_ptsig->Integral());
        //h_ptbg->Scale(1/h_ptbg->Integral());
        float hMax = (h_sig[i]->GetMaximum()>h_bg[i]->GetMaximum()) ? h_sig[i]->GetMaximum() :h_bg[i]->GetMaximum() ;
        h_sig[i]->SetMaximum(hMax*1.1);
        h_sig[i]->SetLineColor(4);
        h_sig[i]->Draw("hist");
        h_bg[i]->Draw("histsame");
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
        /*
        c2->cd();
        h_punziF->Draw("hist");
        c1->cd();
        //h_effiF->Draw("hist");
        TGraph *tg_ptsig = new TGraph(h_sig[i]);
        TGraph *tg_ptbg = new TGraph(h_bg[i]);
        TGraph *tg_effiF = new TGraph(nBin,&info[0][0],&info[1][0]);
        TGraph *tg_effiR = new TGraph(nBin,&info[0][0],&info[2][0]);
        TGraph *tg_punziF = new TGraph(nBin,&info[0][0],&info[3][0]);
        TGraph *tg_punziR = new TGraph(nBin,&info[0][0],&info[4][0]);
        tg_ptsig->Draw("AB");
        */
        Float_t rightmax = 1.1*h_effiF->GetMaximum();
        Float_t scale    = gPad->GetUymax()/rightmax;
        h_effiF->SetLineColor(kRed);
        h_effiF->SetLineWidth(3);
        h_effiF->Scale(scale);
        h_effiF->Draw("histsame");
        //draw an axis on the right side
        TGaxis* axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),0,rightmax,510,"+L");
        axis->SetLineColor(kRed);
        axis->SetLabelColor(kRed);
        axis->Draw();
        
        int lineColor2 = 91;
        Float_t rightmax2 = 1.1*h_punziF->GetMaximum();
        Float_t scale2    = gPad->GetUymax()/rightmax2;
        h_punziF->SetLineColor(lineColor2);
        h_punziF->SetLineWidth(1);
        h_punziF->Scale(scale2);
        h_punziF->Draw("histsame");
        //draw an axis on the right side
        TGaxis* axis2 = new TGaxis(gPad->GetUxmax()*0.9,gPad->GetUymin(),gPad->GetUxmax()*0.9,gPad->GetUymax(),0,rightmax,510,"+L");
        axis2->SetLineColor(lineColor2);
        axis2->SetLabelColor(lineColor2);
        axis2->Draw();
        h_bg[i]->Draw("histsame");
        h_sig[i]->Draw("histsame");
        cout << "i = " << i << "\thName = " << hname[i] << endl;
        c1->Print(pdfName.Data());
    }
    c1->Print((pdfName+"]").Data());
    /*
    vector <double> window = getCut(h_ptsig,h_ptbg);
    cout << "integral:" << h_ptsig->Integral() << endl;
    cout << "window range: " << window[0] << " ~ " << window[1] << endl;
    */
}
