#include <iostream>
#include <fstream>
using namespace std;
void plotEffi() {
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    //c1->SetLogy();
    string mzplist[8] = {"600","800","1000","1200","1400","1700","2000","2500"};
    string njet[4] = {"4","3","3","2"};
    string effi[4][8];
    for (int i=0;i<8;i++) {
        //bool twice=false;
        for (int j=0,k=0;j<3;j++) {
            //if (j==1&&i==7) continue;
            string filename = Form("effi_Zpmass%s_A0mass300_%djets.txt",mzplist[i].data(),4-j);
            ifstream file(filename.data(), ios::in);        //開啟檔案為輸出狀態
            string str;
            while (getline(file,str)) {
                effi[k][i] = str;
                k++;
            }
            file.close();      
        }
    }
    
    vector <TH1F*> h_pass(4);
    vector <TH1F*> h_all(4);
    vector <TGraphAsymmErrors*> gr(4);
    for (int i=0;i<4;i++) {
        h_pass[i] = new TH1F(Form("h_pass_%d",i),"h_pass",10,500,2500);
        h_all[i] = new TH1F(Form("h_all_%d",i),"h_all",10,500,2500);
        gr[i] = new TGraphAsymmErrors();
    }
    TMultiGraph *mg = new TMultiGraph();
    int xbin[8] = {1,2,3,4,5,6,8,10};
    int nPass[8], nFail[8];

    for (int j=0;j<4;j++) {
        for (int i=0;i<8;i++ ) { 
            //if (j==1&&i==7) continue;
            nPass[i] = (int)(strtof((effi[j][i]).c_str(),0)*10000);
            nFail[i] = 10000 - nPass[i];
            h_pass[j]->SetBinContent(xbin[i],nPass[i]);
            h_all[j]->SetBinContent(xbin[i],10000);
        }
        gr[j]->BayesDivide(h_pass[j],h_all[j]);
        gr[j]->SetMarkerStyle(20);
        gr[j]->SetLineWidth(1);
        gr[j]->SetLineColor(j+2);
        gr[j]->SetMarkerColor(j+2);
        if (j==1) gr[j]->SetTitle("3 jets h1a02");
        else if (j==2) gr[j]->SetTitle("3 jets h2a01");
        else gr[j]->SetTitle(Form("%s jets",njet[j].c_str()));
        gr[j]->Draw("ALP");
        mg->Add(gr[j]);
    }
    mg->Draw("ALP");
    mg->SetTitle("efficiency with different analysis");
    mg->GetXaxis()->SetTitle("MZp");
    mg->GetYaxis()->SetTitle("Efficiency");
    c1->BuildLegend(0.15,0.8,0.35,0.9,"","LP");
    c1->Print("effi.pdf");
}
