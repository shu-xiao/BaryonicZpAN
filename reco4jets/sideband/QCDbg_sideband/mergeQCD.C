#include "iostream"
#include "TFile.h"
#include "TKey.h"
#include <TInterpreter.h>
#include <TCanvas.h>
#include <TString.h>
#include "setNCUStyle.C"

#define L2016 35.9*1000 //35.9 fb^-1
//#define L2016 1
//#define drop false

using namespace std;
float getL(int nEvent, float xs) {
    return nEvent/xs;
}
int getNhist(TFile* f) {
    TIter keyList(f->GetListOfKeys()); 
    TKey *key;
    int j = 0;
    while ((key = (TKey*)keyList())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1")) continue;
        j++;
    }
    return j;
}
void mergeQCD(bool drop=true) {
    setNCUStyle(true);

    //bool drop = true; //drop out QCD_HT50to100    
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    // cross-section unit: pb
    const float xsHTbeam[9] = {246400000,27990000,1712000,347700,32100,6831,1207,119.9,25.24};
    string HT_list[9] = {"QCD_HT50to100","QCD_HT100to200","QCD_HT200to300","QCD_HT300to500","QCD_HT500to700","QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf"};
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
    
    const int nHist = getNhist(file[0]); 
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
    TString outputName = (drop)? "QCDbg":"QCDbg_whole";
    c1->Print((outputName+".pdf[").Data());
    
    TLegend *legend = new TLegend(0.78,0.6,0.94,0.83);
    //legend->AddEntry((TObject*)0,"normalized to 2016 L","");
    for (int i=0;i<nHist;i++) { // loop of hist
        TH1F *h_tem[9];
        c1->Clear();
        bool init = true;
        for (int j=0;j<9;j++) { // loop of HT
            if (drop&&j==0) continue; 
            // set base hist
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
            // set the others
            else {
                hmerge[i]->Add(th1f[j][i],L2016/getL(nEvent[j],xsHTbeam[j]));
                h_tem[j] = (TH1F*)hmerge[i]->Clone(hmerge[i]->GetName());
                h_tem[j]->SetLineColor(-j*6+99);
            }
            if (i==HTindex) {
                //h_HTmerge->Add(th1f[j][i],L2016/getL(nEvent[j],xsHTbeam[j]));
            }
            if (i==0) legend->AddEntry(h_tem[j],HT_list[j].data(),"lf");
        }
        int hnum = (drop)? 8:9; 
        // setting histogram
        TString hName = h_tem[8]->GetName();
        static float ymax; 
        ymax = h_tem[8]->GetMaximum()*1.6;
        c1->SetLogy(0);
        if (i==0) h_tem[8]->GetYaxis()->SetTitle("N_{events}");
        else if (hName.Contains("mass")) h_tem[8]->GetXaxis()->SetTitle("M_{jj}^{QCD bg} (GeV)");
        else if (hName.Contains("Pt")) h_tem[8]->GetXaxis()->SetTitle("Pt_{jj}^{QCD bg} (GeV)");
        else if (hName.Contains("Com")) {
            h_tem[8]->GetXaxis()->SetTitle("N_{allCombinations}");
            c1->SetLogy();
            ymax = h_tem[8]->GetMaximum()*100;
        }
        h_tem[8]->GetYaxis()->SetTitle("A.U.");
        h_tem[8]->GetYaxis()->SetTitleOffset(0.5);
        h_tem[8]->GetYaxis()->SetLabelOffset(999);
        h_tem[8]->GetYaxis()->SetLabelSize(0);
        c1->SetLeftMargin(0.1);
        
        h_tem[8]->Draw("hist");
        h_tem[8]->SetMaximum(ymax);
        for (int j=0;j<hnum;j++) h_tem[8-j]->Draw("histsame");
        legend->Draw();
        c1->SetBottomMargin(0.15);
        c1->SetTopMargin(0.15);
        c1->Update();
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
