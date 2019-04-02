#include <iostream>
#include "../setNCUStyle.C"
using namespace std;
void plotHm() {
    setNCUStyle(1);
    gStyle->SetOptTitle(0);
    //vector <int> zplist = {600,800,1000,1200,1400,1700,2000,2500};
    vector <int> zplist = {600,1000,1400,2000};
    vector <TH1F*> histlist;
    TFile *f[10];
    for (int i=0;i<zplist.size();i++) {
        f[i] = TFile::Open(Form("MZP%d_MA0300.root",zplist[i]));
        histlist.push_back((TH1F*)f[i]->Get("h_HmatchM"));
    }

    TCanvas *c1 = new TCanvas("c1","c1",800,700);
    TLegend *leg = new TLegend(0.65,0.7,0.93,0.9);
    TLatex tex;
    tex.SetTextSize(0.04);
    for (int i=0;i<zplist.size();i++) {
        histlist[i]->SetLineColor(99-7*i);
        histlist[i]->Scale(1/histlist[i]->GetEntries());
        leg->AddEntry(histlist[i],Form("M_{Z'} = %d GeV",zplist[i] ),"al");
    }
    histlist[0]->SetMaximum(0.03);
    histlist[0]->GetXaxis()->SetTitle("M_{h} (GeV)");
    histlist[0]->Draw("hist");
    histlist[0]->GetYaxis()->SetLabelSize(0.04);
    histlist[0]->GetYaxis()->SetTitle("normalized to 1");
    for (int i=1;i<zplist.size();i++) {
        histlist[i]->Draw("histsame");
    }
    leg->Draw();
    tex.DrawLatexNDC(0.7,0.65,"M_{A_{0}} = 300 GeV");
    c1->SetLeftMargin(0.2);
    c1->Print("higgsmass.pdf");
}
