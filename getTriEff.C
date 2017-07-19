#include "iostream"
#include <fstream>
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
using namespace std;
float getTriEff_base(string FName) {
    fstream file;
    file.open(FName.data(),ios::in);
    //string line;
    char line[10];
    //file.read(line.data(),8)
    file >> line;
    file.close();
    return stof(line);
}
void getTriEff() {
    TCanvas *c1 = new TCanvas("c1","c1", 1600, 900);
    c1->SetGrid();
    TH1F *h_pass = new TH1F("h_pass","h_pass",20,0,2100);
    TH1F *h_all = new TH1F("h_all","h_all",20,0,2100);
    TGraphAsymmErrors *gr = new TGraphAsymmErrors();
    float triEff[20], mZp[20], ex[20], ey[20];
    int nPass[20], nFail[20]; 
    for (int i=1;i<=20;i++ ) { 
        triEff[i] = getTriEff_base(Form("triEff/triEff_MZp%d.txt",i*100));
        nPass[i] = (int)(triEff[i]*10000);
        nFail[i] = 10000 - nPass[i];
        mZp[i] = i*100;
        h_pass->SetBinContent(i,nPass[i]);
        h_all->SetBinContent(i,10000);
        //h_pass->Fill(i*100,nPass[i]);
        //h_all->Fill(i*100,10000);
    }
    //for (int j=0;j<20;j++) cout << triEff[j] << endl;
    //gr = new TGraph(20,mZp,triEff);
    //gr = new TGraphErrors(20,mZp,triEff,ex,ey);
    gr->BayesDivide(h_pass,h_all);
    /*
    gr->SetLineColor(1);
    gr->SetLineWidth(2);
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerStyle(22);
    gr->SetMarkerSize(3);
    gr->SetTitle("trigger efficiency with MZp");
    gr->GetYaxis()->SetTitle("efficiency");
    gr->GetXaxis()->SetTitle("MZp (GeV)");
    */
    gr->Draw("AP");
    c1->SaveAs("~/triEff.png");
    c1->Print("triEff.pdf");
}

