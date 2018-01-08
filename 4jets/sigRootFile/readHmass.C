#include <iostream>
#include "TFile.h"
#include "TCanvas.h"

using namespace std;
void norm(TH1F* h1) {
    float total = h1->Integral();
    h1->Scale(1/total);
}
void readHmass() {
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    int ZpMass[8] = {600,800,1000,1200,1400,1700,2000,2500};
    TH1F *hMass[8];
    TH1F *hMass_genMatch[8];
    for (int i=0;i<8;i++) {
        TFile *f = TFile::Open(Form("bb2HDM_MZp%d_MA0300.root",ZpMass[i]));
        TFile *f1 = TFile::Open(Form("bb2HDM_genMatch_MZp%d_MA0300.root",ZpMass[i]));
        hMass[i] =  (TH1F*)f->Get("h_higgsM_effi");
        hMass_genMatch[i] =  (TH1F*)f1->Get("h_higgsM_ori");
        norm(hMass[i]);
        norm(hMass_genMatch[i]);
    }
    hMass[3]->SetTitle("h_higgsM_applyingSelection");
    hMass[3]->SetLineColor(99-30);
    hMass[3]->Draw("hist");
    for (int i=2;i>=0;i--) {
        hMass[i]->SetLineColor(99-i*10);
        //hMass[i]->DrawNormalized("histsame");
        hMass[i]->Draw("histsame");
    }
    
    TLegend* legend = new TLegend(0.2,0.7,0.4,0.9);
    for (int i=0;i<4;i++) legend->AddEntry(hMass[i],Form("MZp %d GeV",ZpMass[i]),"l");
    legend->Draw();
    c1->Print("hMass_MZp600to1200.pdf(");
    
    hMass_genMatch[3]->SetTitle("h_higgsM_genMatch");
    hMass_genMatch[3]->SetLineColor(99-30);
    hMass_genMatch[3]->Draw("hist");
    for (int i=2;i>=0;i--) {
        hMass_genMatch[i]->SetLineColor(99-i*10);
        //hMass[i]->DrawNormalized("histsame");
        hMass_genMatch[i]->Draw("histsame");
    }
    TLegend* legend2 = new TLegend(0.2,0.7,0.4,0.9);
    for (int i=0;i<4;i++) legend2->AddEntry(hMass_genMatch[i],Form("MZp %d GeV",ZpMass[i]),"l");
    legend2->Draw();
    c1->Print("hMass_MZp600to1200.pdf)");
}
