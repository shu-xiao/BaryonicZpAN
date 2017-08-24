#include "iostream"
#include "TFile.h"
#include "TKey.h"
#include <TInterpreter.h>
#include <TCanvas.h>
#include <TString.h>
#define L2016 35.9*1000 //35.9 fb^-1
//#define L2016 1
using namespace std;
float getL(int nEvent, float xs) {
    return nEvent/xs;
}
void mergeQCD() {
    
    bool drop = true; //drop out QCD_HT50to100    
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    // cross-section unit: pb
    float xsHTbeam[9] = {246400000,27990000,1712000,347700,32100,6831,1207,119.9,25.24};
    TFile* file[9];
    file[0] = TFile::Open("QCD_HT50to100.root");
    file[1] = TFile::Open("QCD_HT100to200.root");
    file[2] = TFile::Open("QCD_HT200to300.root");
    file[3] = TFile::Open("QCD_HT300to500.root");
    file[4] = TFile::Open("QCD_HT500to700.root");
    file[5] = TFile::Open("QCD_HT700to1000.root");
    file[6] = TFile::Open("QCD_HT1000to1500.root");
    file[7] = TFile::Open("QCD_HT1500to2000.root");
    file[8] = TFile::Open("QCD_HT2000toInf.root");
    
    const int nHist = 25; 
    TH1F* hmerge[nHist];
    TH1F* th1f[9][nHist];
    TH1F* h_HTmerge = new TH1F();
    h_HTmerge->Sumw2();
    vector <int> nTii, nEvent;
    int HTindex = -1; 
    for (int i=0;i<9;i++) {
        TIter keyList(file[i]->GetListOfKeys()); 
        TKey *key;
        int j = 0;
        while ((key = (TKey*)keyList())) {
            TClass *cl = gROOT->GetClass(key->GetClassName());
            if (!cl->InheritsFrom("TH1")) continue;
            th1f[i][j] = (TH1F*)key->ReadObj();
            TString hName = th1f[i][j]->GetName();
            if ( hName.Contains("h_allEvent")) {
                nEvent.push_back(th1f[i][j]->GetEntries());
            }
            if (hName.Contains("h_HT")) HTindex = j;
            j++;
        }
    }
    cout << HTindex << endl;
    TString outputName = (drop)? "QCDbg":"QCDbg_whole";
    c1->Print((outputName+".pdf[").Data());
    for (int i=0;i<nHist;i++) {
        TH1F *h_tem[9];
        c1->Clear();
        bool init = true;
        for (int j=0;j<9;j++) {
            if (drop&&j==0) continue; 
            if (init) {
                hmerge[i] = (TH1F*)th1f[j][i]->Clone(th1f[j][i]->GetName());
                hmerge[i]->Sumw2();
                hmerge[i]->Scale(L2016/getL(nEvent[j],xsHTbeam[j]));
                h_tem[j] = (TH1F*)th1f[j][i]->Clone(th1f[j][i]->GetName());
                h_tem[j]->Sumw2();
                h_tem[j]->Scale(L2016/getL(nEvent[j],xsHTbeam[j]));
                h_tem[j]->SetLineColor(99);
                init = false;
            }
            else {
                hmerge[i]->Add(th1f[j][i],L2016/getL(nEvent[j],xsHTbeam[j]));
                h_tem[j] = (TH1F*)hmerge[i]->Clone(hmerge[i]->GetName());
                h_tem[j]->SetLineColor(-j*6+99);
            }
            if (i==HTindex) {
                //h_HTmerge->Add(th1f[j][i],L2016/getL(nEvent[j],xsHTbeam[j]));
            }
        }
        int hnum = (drop)? 8:9; 
        h_tem[8]->Draw("hist");
        if (i==HTindex||i==0) c1->SetLogy();
        else c1->SetLogy(0);
        for (int j=0;j<hnum;j++) h_tem[8-j]->Draw("histsame");
        c1->Print((outputName+".pdf").Data());
    }
    //h_HTmerge->Draw("hist");
    //c1->Print((outputName+".pdf").Data());
    c1->Print((outputName+".pdf]").Data());
    //for (int i=0;i<nTi.size();i++) cout << nTi[i] << " ";
    /*
    for (int i=0;i<nHist;i++) {
        hmerge[i]->Draw("hist");
        th1f[0][i]->SetLineColor(4);
        th1f[0][i]->Scale(L2016/getL(nEvent[0],xsHTbeam[0]));
        th1f[0][i]->Draw("histsame");
        c1->Print("QCDbg.pdf");
    }
    */
    
    TFile* output = new TFile((outputName+".root").Data(),"recreate");
    for (int i=0;i<nHist;i++) hmerge[i]->Write();
    //hmerge[0]->Write();
    output->Close();
}
