#include <iostream>

void fitWidth() {
    TFile *f = TFile::Open("fitwidth.root");
    TH1F *h_hM = (TH1F*)f->Get("h_higgsM");
    TH1F *h_a0M = (TH1F*)f->Get("h_a0M");
    TH1F *h_zpM = (TH1F*)f->Get("h_ZpM");
    TF1 *fit = new TF1("fit","crystalball");
    float p1, p2, p3, p4, p5;
    fit->SetParameters(p1,p2,p3,p4,p5);
    h_hM->Fit(fit);

    
}
