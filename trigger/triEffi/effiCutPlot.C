#include <iostream>
#include <fstream>

using namespace std;
int* nEventList (string fileName) {
    int* data;
    const int size = 6;
    ifstream inputFile(fileName.data());
    data = new int[size];  
    for(int i=0;i<size;i++)
    {
        inputFile >> data[i];
    }
    return data;
}
void effiCutPlot() {
    gStyle->SetOptStat(0);
    
    TCanvas *c1 = new TCanvas("c1","c1",1100,600);
    gPad->SetGrid();
    const int nCut = 4;
    int zpMass[8] = {600,800,1000,1200,1400,1700,2000,2500};
    float zprange[9] = {500,700,900,1100,1300,1500,1800,2200,2800};
    TH1F* h_passTri[2];
    for (int i=0;i<2;i++) h_passTri[i] = new TH1F(Form("h_passTri%d",i),Form("h_passTri%d",i),8,zprange);
    TH1F* h_pass = new TH1F("h_pass","h_pass",8,zprange);
    TH1F* h_passHT = new TH1F("h_passHT","h_passtrigger||HT250trigger",8,zprange);
    TH1F* h_passHT800 = new TH1F("h_passHT","h_passtrigger||HT800trigger",8,zprange);
    TH1F* h_all = new TH1F("h_all","h_all",8,zprange);
    for (int i=0;i<2;i++) h_passTri[i]->Sumw2();
    h_pass->Sumw2();
    h_passHT->Sumw2();
    h_passHT800->Sumw2();
    h_all->Sumw2();
    int xbin[8] = {1,2,3,4,5,6,8,10};
    for (int n=0;n<8;n++) {
        string fname = Form("triEffi_MZp%d_MA0300.txt",zpMass[n]); 
        int* nEvent = nEventList(fname);
        for (int nE=0;nE<nEvent[0];nE++) h_all->Fill(zpMass[n]);
        for (int i=0;i<2;i++) for (int nE=0;nE<nEvent[i+1];nE++)  h_passTri[i]->Fill(zpMass[n]);
        for (int nE=0;nE<nEvent[3];nE++) h_pass->Fill(zpMass[n]);
        for (int nE=0;nE<nEvent[4];nE++) h_passHT->Fill(zpMass[n]);
        for (int nE=0;nE<nEvent[5];nE++) h_passHT800->Fill(zpMass[n]);
        //h_all->Fill(zpMass[n],nEvent[0]);
        //for (int i=0;i<2;i++) h_passTri[i]->Fill(zpMass[n],nEvent[i+1]);
        //h_pass->Fill(zpMass[n],nEvent[3]);
        for (int i=0;i<4;i++) {
            //cout << *nEvent << endl;
            nEvent++;
        }
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
    TGraphAsymmErrors* gr = new TGraphAsymmErrors(h_pass,h_all,"cl=0.683 b(1,1) mode");
    TGraphAsymmErrors* grHT = new TGraphAsymmErrors(h_passHT,h_all,"cl=0.683 b(1,1) mode");
    TGraphAsymmErrors* grHT800 = new TGraphAsymmErrors(h_passHT800,h_all,"cl=0.683 b(1,1) mode");
    TGraphAsymmErrors* gr_Tri[2];
    for (int i=0;i<2;i++) gr_Tri[i] = new TGraphAsymmErrors(h_passTri[i],h_all,"cl=0.683 b(1,1) mode");
    gr->GetYaxis()->SetRangeUser(0., 1.);
    for (int i=0;i<2;i++) gr_Tri[i]->GetYaxis()->SetRangeUser(0., 1.);
    gr->SetMarkerStyle(20);
    gr->SetLineWidth(2);
    gr->GetXaxis()->SetTitle("MZp (GeV)");
    gr->GetYaxis()->SetTitle("Efficiency");
    gr->SetLineColor(1);
    gr->SetMarkerColor(1);
    gr->SetTitle("Trigger Efficiency");
    TLatex latex;
    latex.SetTextSize(0.025);
    latex.SetTextAlign(23);
    
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
    gr->Draw("ALP");
    c1->Update();
    c1->Print("effiCutPlot.pdf(");
    grHT->Draw("ALP");
    c1->Print("effiCutPlot.pdf");
    grHT800->Draw("ALP");
    c1->Print("effiCutPlot.pdf");
    gr_Tri[0]->Draw("ALP");
    c1->Print("effiCutPlot.pdf");
    gr_Tri[1]->Draw("ALP");
    c1->Print("effiCutPlot.pdf)");
}
