#include <iostream>
#include <fstream>
using namespace std;
void plotEffi() {
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    string mzplist[8] = {"600","800","1000","1200","1400","1700","2000","2500"};
    string njet[3] = {"4","3","2"};
    string effi[3][8];
    string str;
    for (int i=0;i<8;i++) {
        for (int j=0;j<3;j++) {
            if (j==1&&i==7) continue;
            string filename = Form("effi_Zpmass%s_A0mass300_%sjets.txt",mzplist[i].data(),njet[j].data());
            ifstream file(filename.data(), ios::in);        //開啟檔案為輸出狀態
            while (getline(file,str)) {
                effi[j][i] = str;
                break;
                //cout << str << endl;
            }
            file.close();      
        }
    }

    vector <TH1F*> h_pass(3);
    vector <TH1F*> h_all(3);
    vector <TGraphAsymmErrors*> gr(3);
    for (int i=0;i<3;i++) {
        h_pass[i] = new TH1F(Form("h_pass_%d",i),"h_pass",10,500,2500);
        h_all[i] = new TH1F(Form("h_all_%d",i),"h_all",10,500,2500);
        gr[i] = new TGraphAsymmErrors();
    }
    TMultiGraph *mg = new TMultiGraph();
    int xbin[8] = {1,2,3,4,5,6,8,10};
    int nPass[8], nFail[8];
    for (int j=0;j<3;j++) {
        for (int i=0;i<8;i++ ) { 
            if (j==1&&i==7) continue;
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
        gr[j]->SetTitle(Form("%d jets",j+2));
        gr[j]->Draw("ALP");
        mg->Add(gr[j]);
    }
    mg->Draw("ALP");
    mg->SetTitle("efficiency with different analysis");
    mg->GetXaxis()->SetTitle("MZp");
    mg->GetYaxis()->SetTitle("Efficiency");
    c1->BuildLegend(0.7,0.8,0.9,0.9,"","LP");
    c1->Print("effi.pdf");
}
