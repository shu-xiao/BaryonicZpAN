#include "iostream"
#include "TFile.h"
#include "TKey.h"
#include <TInterpreter.h>
#include <TCanvas.h>
#include <TString.h>

#define L2016 35.9*1000 //35.9 fb^-1
using namespace std;
float getL(int nEvent, float xs) {
    return nEvent/xs;
}
void mergeQCD() {
    
    
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
    
    const int nHist = 19; 
    TH1F *hmerge[nHist];
    TH1F* th1f[9][nHist];
    vector <int> nTii, nEvent;
     
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
            j++;
        }
    }
    for (int i=0;i<nHist;i++) {
        for (int j=0;j<9;j++) {
            if (j==0) {
                hmerge[i] = (TH1F*)th1f[j][i]->Clone(th1f[j][i]->GetName());
                hmerge[i]->Sumw2(L2016/getL(nEvent[j],xsHTbeam[j]));
            }
            else hmerge[i]->Add(th1f[j][i],L2016/getL(nEvent[j],xsHTbeam[j]));
        }
    }
    //for (int i=0;i<nTi.size();i++) cout << nTi[i] << " ";
    c1->Print("QCDbg.pdf[");
    for (int i=0;i<nHist;i++) {
        hmerge[i]->Draw("hist");
        c1->Print("QCDbg.pdf");
    }
    c1->Print("QCDbg.pdf]");
    TFile* output = new TFile("QCDbg.root","recreate");
    for (int i=0;i<nHist;i++) hmerge[i]->Write();
    //hmerge[0]->Write();
    output->Close();
    cout << "finish !" << endl;
}
