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
    int zpMass[8] = {600,800,1000,1200,1400,1700,2000,2500};
    TH1F* effiPlot[7];
    for (int i=0;i<7;i++) effiPlot[i] = new TH1F(Form("effiPlot%d",i),"Efficiency of Four Jets Category Varies with MZp and Selections shown with Colors",10,0,10);
    vector <int*> eventList;
    eventList.reserve(8);
    for (int n=0;n<8;n++) {
        string fname = Form("effi_Zpmass%d_A0mass300_4jets.txt",zpMass[n]); 
        int* nEvent = nEventList(fname);
        /*
        for (int i=0;i<6;i++) {
            cout << *nEvent << endl;
            nEvent++;
        }
        */
        eventList.push_back(nEvent);
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
    for (int i=0;i<8;i++) { // different mass sample
        int *nEvent = eventList.at(i);
        //effiPlot[i]->SetLineWdith(2);
        
        for (int j=0;j<6;j++) { // different cut
             effiPlot[j]->Fill(Form("%d",zpMass[i]),*nEvent);
             nEvent++;
        }
        
    }
    effiPlot[0]->Draw("hist");
    effiPlot[0]->GetXaxis()->SetTitle("MZp (GeV)");
    effiPlot[0]->GetYaxis()->SetTitle("number of events");
    effiPlot[0]->GetYaxis()->SetTitleOffset(1.2);
    int color;
    for (int i=0;i<6;i++) 
    {
        effiPlot[i]->Fill("",0);
        effiPlot[i]->Fill("color tickets",10000*(6-i)/6);
        color = 99-7*i;
        //effiPlot[i]->SetLineColor(1);
        effiPlot[i]->SetLineColor(color);
        effiPlot[i]->SetFillColor(color);
        effiPlot[i]->Draw("histsame");

    }
    TLatex latex;
    latex.SetTextSize(0.025);
    latex.SetTextAlign(23);
    float y = 10000/12;
    latex.DrawLatex(9.5,11*y,"All Event");
    latex.DrawLatex(9.5,9*y,"njets>=4");
    latex.DrawLatex(9.5,7*y,"110<M_{h}<140");
    latex.DrawLatex(9.5,5*y,"250<M_{A0}<350");
    latex.DrawLatex(9.5,3*y,"#splitline{|M_{4jets}-M_{Zp}|}{<100}");
    latex.DrawLatex(9.5,1*y,"N_{candidates}=1");
    //latex.DrawLatex(8.5,1*y,"#splitline{left event}{cut =>}");
    //c1->Update();
    c1->Print("effiCutPlot.pdf");
}
