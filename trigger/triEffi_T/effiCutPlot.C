#include <iostream>
#include <fstream>
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TH1.h"
#include "TCanvas.h"

using namespace std;
int* nEventList (string fileName) {
    int* data;
    const int size = 7;
    ifstream inputFile(fileName.data());
    data = new int[size];  
    for(int i=0;i<size;i++)
    {
        inputFile >> data[i];
    }
    return data;
}
void getEventList(string fileName, int arr1[],int arr2[]) {
    const int size = 8;
    ifstream inputFile(fileName.data());
    for(int i=0;i<size;i++)
    {
        inputFile >> arr1[i] >> arr2[i];
    }
    
}
void effiCutPlot() {
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(1);
    gStyle->SetMarkerStyle(20);
    
    int tagColor = 2;
    int antiTagColor = 4;
    int Lwid = 1;
    
    TCanvas *c1 = new TCanvas("c1","c1",1100,600);
    gPad->SetGrid();
    const int nCut = 4;
    int zpMass[8] = {600,800,1000,1200,1400,1700,2000,2500};
    float zprange[9] = {500,700,900,1100,1300,1500,1800,2200,2800};
    TH1F* h_passTri[2];
    TH1F* h_passTri_anti[2];
    for (int i=0;i<2;i++) h_passTri[i] = new TH1F(Form("h_passTri%d",i),"h_Tag",8,zprange);
    for (int i=0;i<2;i++) h_passTri_anti[i] = new TH1F(Form("h_passTri%d_anti",i),"h_antiTag",8,zprange);
    TH1F* h_pass = new TH1F("h_tag","h_Tag",8,zprange);
    TH1F* h_pass_anti = new TH1F("h_antitag","h_antiTag",8,zprange);
    TH1F* h_all = new TH1F("h_all","h_all",8,zprange);
    TH1F* h_all_anti = new TH1F("h_all_anti","h_all_anti",8,zprange);
    for (int i=0;i<2;i++) h_passTri[i]->Sumw2();
    h_pass->Sumw2();
    h_all->Sumw2();
    int xbin[8] = {1,2,3,4,5,6,8,10};
    for (int n=0;n<8;n++) {
        string fname = Form("triEffi_MZp%d_MA0300.txt",zpMass[n]);
        int nEvent[8];
        int nEvent2[8];
        for (int i=0;i<8;i++) nEvent[i]=0;
        for (int i=0;i<8;i++) nEvent2[i]=0;
        //int* nEvent = nEventList(fname);
        getEventList(fname,nEvent,nEvent2);
        for (int nE=0;nE<nEvent[0];nE++) h_all->Fill(zpMass[n]);
        for (int i=0;i<2;i++) for (int nE=0;nE<nEvent[i+2];nE++)  h_passTri[i]->Fill(zpMass[n]);
        for (int nE=0;nE<nEvent[4];nE++) h_pass->Fill(zpMass[n]);
        
        for (int nE=0;nE<nEvent2[0];nE++) h_all_anti->Fill(zpMass[n]);
        for (int i=0;i<2;i++) for (int nE=0;nE<nEvent2[i+2];nE++)  h_passTri_anti[i]->Fill(zpMass[n]);
        for (int nE=0;nE<nEvent2[4];nE++) h_pass_anti->Fill(zpMass[n]);
    }
    /*
    for (int i=0;i<8;i++) {
        for (int j=0;j<6;j++) {
            cout << *eventList.at(i) << "\t";
            eventList.at(i)++;
        }
        cout << endl;
    }
    */
    TMultiGraph* mg1 = new TMultiGraph("mg1","Trigger Efficiency (tr1 or tr2)");
    TMultiGraph* mg2 = new TMultiGraph("mg2","tr1 (DoubleJet90_Double30_TripleBTagCSV_p087)");
    TMultiGraph* mg3 = new TMultiGraph("mg3","tr2 (QuadJet45_TripleBTagCSV_p087)");
    TMultiGraph* mglist[2] = {mg2,mg3};
    TGraphAsymmErrors* gr = new TGraphAsymmErrors(h_pass,h_all,"cl=0.683 b(1,1) mode");
    TGraphAsymmErrors* gr2 = new TGraphAsymmErrors(h_pass_anti,h_all_anti,"cl=0.683 b(1,1) mode");
    
    TGraphAsymmErrors* gr_Tri[2];
    TGraphAsymmErrors* gr_Tri_anti[2];
    for (int i=0;i<2;i++) gr_Tri[i] = new TGraphAsymmErrors(h_passTri[i],h_all,"cl=0.683 b(1,1) mode");
    //gr->GetYaxis()->SetRangeUser(0., 1.);
    //for (int i=0;i<2;i++) gr_Tri[i]->GetYaxis()->SetRangeUser(0., 1.);
    for (int i=0;i<2;i++) gr_Tri_anti[i] = new TGraphAsymmErrors(h_passTri_anti[i],h_all_anti,"cl=0.683 b(1,1) mode");
    //gr2->GetYaxis()->SetRangeUser(0., 1.);
    //for (int i=0;i<2;i++) gr_Tri_anti[i]->GetYaxis()->SetRangeUser(0., 1.);
    
    gr->SetLineColor(tagColor);
    gr->SetMarkerColor(tagColor);
    gr2->SetLineColor(antiTagColor);
    gr2->SetMarkerColor(antiTagColor);
    gr->SetLineWidth(Lwid);
    gr2->SetLineWidth(Lwid);
    mg1->Add(gr);
    mg1->Add(gr2);
    mg1->Draw("APL");
    c1->Update();
    mg1->GetYaxis()->SetRangeUser(0., 1.);
    mg1->GetXaxis()->SetTitle("MZp (GeV)");
    mg1->GetYaxis()->SetTitle("Efficiency");
    TLatex latex;
    latex.SetTextSize(0.025);
    latex.SetTextAlign(23);
    
    for (int i=0;i<2;i++) {
        gr_Tri[i]->SetLineWidth(Lwid);
        gr_Tri_anti[i]->SetLineWidth(Lwid);
        mglist[i]->Add(gr_Tri[i]);
        mglist[i]->Add(gr_Tri_anti[i]);
        mglist[i]->Draw("APL");
        gr_Tri[i]->SetLineColor(tagColor);
        gr_Tri[i]->SetMarkerColor(tagColor);
        gr_Tri_anti[i]->SetLineColor(antiTagColor);
        gr_Tri_anti[i]->SetMarkerColor(antiTagColor);
        c1->Update();
        mglist[i]->GetYaxis()->SetRangeUser(0., 1.);
        mglist[i]->GetXaxis()->SetTitle("MZp (GeV)");
        mglist[i]->GetYaxis()->SetTitle("Efficiency");
    }
    
    /*
    string label[nCut] = {"All Events","nAK4jet>=4","nGoodJet>=4","#chi^{2}<100"};
    for (int i=0;i<nCut;i++) {
        latex.DrawLatex(9.5,(nCut*2-1-2*i)*y,label[i].data());
    }
    latex.DrawLatex(9.5,11*y,"All Event");
    latex.DrawLatex(9.5,9*y,"njets>=4");
    latex.DrawLatex(9.5,7*y,"110<M_{h}<140");
    latex.DrawLatex(9.5,5*y,"250<M_{A0}<350");
    latex.DrawLatex(9.5,3*y,"#splitline{|M_{4jets}-M_{Zp}|}{<100}");
    latex.DrawLatex(9.5,1*y,"N_{candidates}=1");
    */
    //latex.DrawLatex(8.5,1*y,"#splitline{left event}{cut =>}");
    //c1->Update();
    mg1->Draw("ALP");
    c1->BuildLegend();
    //c1->Update();
    c1->Print("effiCutPlot.pdf(");
    mglist[0]->Draw("ALP");
    c1->BuildLegend();
    c1->Print("effiCutPlot.pdf");
    mglist[1]->Draw("ALP");
    c1->BuildLegend();
    c1->Print("effiCutPlot.pdf)");
}
