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
void mergeQCD(bool drop=true) {
    setNCUStyle(true);

    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    // cross-section unit: pb
    const float xsHTbeam[9] = {246400000,27990000,1712000,347700,32100,6831,1207,119.9,25.24};
    string HT_list[9] = {"QCD_TMVA_HT50to100","QCD_TMVA_HT100to200","QCD_TMVA_HT200to300","QCD_TMVA_HT300to500","QCD_TMVA_HT500to700","QCD_TMVA_HT700to1000","QCD_TMVA_HT1000to1500","QCD_TMVA_HT1500to2000","QCD_TMVA_HT2000toInf"};
    TFile* file[9];
    file[0] = TFile::Open("QCD_HT50to100_T2.root");
    file[1] = TFile::Open("QCD_TMVA_HT100to200.root");
    //file[1] = TFile::Open("QCD_HT100to200_T2.root");
    file[2] = TFile::Open("QCD_HT200to300_T2.root");
    file[3] = TFile::Open("QCD_HT300to500_T2.root");
    file[4] = TFile::Open("QCD_HT500to700_T2.root");
    file[5] = TFile::Open("QCD_HT700to1000_T2.root");
    file[6] = TFile::Open("QCD_HT1000to1500_T2.root");
    file[7] = TFile::Open("QCD_HT1500to2000_T2.root");
    file[8] = TFile::Open("QCD_HT2000toInf_T2.root");
    /*
    file[0] = TFile::Open("QCD_TMVA_HT50to100.root");
    file[1] = TFile::Open("QCD_TMVA_HT100to200.root");
    file[2] = TFile::Open("QCD_TMVA_HT200to300.root");
    file[3] = TFile::Open("QCD_TMVA_HT300to500.root");
    file[4] = TFile::Open("QCD_TMVA_HT500to700.root");
    file[5] = TFile::Open("QCD_TMVA_HT700to1000.root");
    file[6] = TFile::Open("QCD_TMVA_HT1000to1500.root");
    file[7] = TFile::Open("QCD_TMVA_HT1500to2000.root");
    file[8] = TFile::Open("QCD_TMVA_HT2000toInf.root");
    */
    // empty 3 4 0 5
    // load TH1F and TH2D from root files 
    vector <TH1F*> hmerge;
    vector <vector<TH1F*>> th1f(9);
    vector <vector<TH2F*>> th2f(9);
    vector <TH2F*> th2f_sum;
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
            // TH2D
            if (cl->InheritsFrom("TH2F")) {
                th2f[i].push_back((TH2F*)key->ReadObj()); 
            }
            if (!cl->InheritsFrom("TH1")||cl->InheritsFrom("TH2F")) continue;
            // TH1F
            TH1F* temth1f = (TH1F*)key->ReadObj();
            TString hName = temth1f->GetName();
            if ( hName.Contains("h_allEvent")) {
                nEvent.push_back(temth1f->GetEntries());
            }
            if (hName.Contains("h_HT")) HTindex = j;
            th1f[i].push_back(temth1f);
            j++;
        }
    }
    // normalize to 2017 luminosity
    TString outputName = (drop)? "QCDbg":"QCDbg_whole";
    c1->Print((outputName+".pdf[").Data());
    
    TLegend *legend = new TLegend(0.78,0.6,0.94,0.83);
    //legend->AddEntry((TObject*)0,"normalized to 2016 L","");
    cout << "merge TH1F" << endl;
    for (int i=0;i<th1f[0].size();i++) { // loop of hist
        TH1F *h_tem[9];
        c1->Clear();
        bool init = true;
        // copy first file
        int iniInd = (drop) ?1:0;
        /*
        if (i==3||i==4||i==5) {
            hmerge.push_back(new TH1F());
            continue;
        }
        */
        hmerge.push_back((TH1F*)th1f[iniInd][i]->Clone(th1f[iniInd][i]->GetName()));
        hmerge[i]->Sumw2();
        if (nEvent[iniInd]) hmerge[i]->Scale(L2016/getL(nEvent[iniInd],xsHTbeam[iniInd]));
        h_tem[iniInd] = (TH1F*)th1f[iniInd][i]->Clone(th1f[iniInd][i]->GetName());
        h_tem[iniInd]->Sumw2();
        if (nEvent[iniInd]) h_tem[iniInd]->Scale(L2016/getL(nEvent[iniInd],xsHTbeam[iniInd]));
        h_tem[iniInd]->SetLineColor(99);
        for (int j=iniInd;j<9;j++) { // loop of HT
            //if (drop&&j==0) continue; 
            // set the others
            if (nEvent[j]) hmerge[i]->Add(th1f[j][i],L2016/getL(nEvent[j],xsHTbeam[j]));
            h_tem[j] = (TH1F*)hmerge[i]->Clone(hmerge[i]->GetName());
            h_tem[j]->SetLineColor((-j+iniInd)*6+99);
            //h_tem[j]->SetLineColor((-j+iniInd)*6+99);
            if (i==HTindex) {
                //h_HTmerge->Add(th1f[j][i],L2016/getL(nEvent[j],xsHTbeam[j]));
            }
            if (i==0) legend->AddEntry(h_tem[j],HT_list[j].data(),"lf");
        }
        // setting histogram
        int hnum = (drop)? 8:9; 
        TString hName = h_tem[8]->GetName();
        static float ymax;
        ymax = 0;
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
        //for (int j=0;j<hnum;j++) if (8-j!=3||8-j!=4||8-j!=5)h_tem[8-j]->Draw("histsame");
        for (int j=0;j<hnum;j++) h_tem[8-j]->Draw("histsame");
        legend->Draw();
        //c1->SetBottomMargin(0.15);
        //c1->SetTopMargin(0.15);
        c1->Update();
        c1->Print((outputName+".pdf").Data());
    }
    //h_HTmerge->Draw("hist");
    
    // setting th2f
    cout << "merge TH2F" << endl;
    gStyle->SetOptTitle(0);
    for (int i=0;i<th2f[0].size();i++) {
        c1->Clear();
        int iniInd = (drop)?1:0;
        TH2F* th2f_tem = (TH2F*)th2f[iniInd][i]->Clone(th2f[iniInd][i]->GetName());
        th2f_tem->Sumw2();
        if (nEvent[iniInd]) th2f_tem->Scale(L2016/getL(nEvent[iniInd],xsHTbeam[iniInd]));
        
        for (int j=1;j<9;j++) {
            if (drop&&j==1) continue;
            if (nEvent[j]) th2f_tem->Add(th2f[j][i],L2016/getL(nEvent[j],xsHTbeam[j]));
        }
        // setting hist
        TString name = th2f_tem->GetName(); 
        if (name.Contains("ha0")) {
            th2f_tem->GetXaxis()->SetTitle("(M_{jj}-M_{h})/#sigma_{h}");
            th2f_tem->GetYaxis()->SetTitle("(M_{jj}-M_{A0})/#sigma_{A0}");
        }
        if (name.Contains("hzp")) {
            th2f_tem->GetXaxis()->SetTitle("(M_{jj}-M_{h})/#sigma_{h}");
            th2f_tem->GetYaxis()->SetTitle("(M_{4j}-M_{Zp})/#sigma_{Zp}");
        }
        if (name.Contains("a0zp")) {
            th2f_tem->GetXaxis()->SetTitle("(M_{jj}-M_{A0})/#sigma_{A0}");
            th2f_tem->GetYaxis()->SetTitle("(M_{4j}-M_{Zp})/#sigma_{Zp}");
        }
        th2f_tem->SetTitleFont(62,"XYZ");
        th2f_sum.push_back(th2f_tem);
        th2f_tem->Draw("colz");
        c1->SetLeftMargin(0.15);
        c1->SetRightMargin(0.15);
        c1->SetBottomMargin(0.15);
        c1->Update();
        c1->Print((outputName+".pdf").Data());
    }
    
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
    
    // write root file
    TFile* output = new TFile((outputName+".root").Data(),"recreate");
    for (int i=0;i<hmerge.size();i++) hmerge[i]->Write();
    for (int i=0;i<th2f_sum.size();i++) th2f_sum[i]->Write();
    output->Close();
}
