#include <iostream>
#include "TFile.h"
#include "TCanvas.h"

using namespace std;
void readHmass() {
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    int ZpMass[8] = {600,800,1000,1200,1400,1700,2000,2500};
    TH1F *hMass[8];
    for (int i=0;i<8;i++) {
        TFile *f = TFile::Open(Form("bb2HDM_MZp%d_MA0300.root",ZpMass[i]));
        hMass[i] =  (TH1F*)f->Get("h_higgsM");
    }
    hMass[3]->Draw("hist");
    hMass[3]->SetLineColor(99-30);
    for (int i=2;i>=0;i--) {
        hMass[i]->Draw("histsame");
        hMass[i]->SetLineColor(99-i*10);
    }
    
    TLegend* legend = new TLegend(0.2,0.7,0.4,0.9);
    for (int i=0;i<4;i++) legend->AddEntry(hMass[i],Form("MZp %d GeV",ZpMass[i]),"l");
    legend->Draw();
    c1->Print("hMass_MZp600to1200.pdf");
}
