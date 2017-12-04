
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1D.h>
#include <TFile.h>
#include "../untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <fstream>
#include <TCanvas.h>
#include "string"
#include "../setNCUStyle.C"
#include "genMatch.C"
using namespace std;
void test(vector<vector<vector<int>>> matchjetlist,int mode=0) {
    int nMatch = 0;
        for (int i=0;i<matchjetlist.size();i++) {
            if (matchjetlist[i].size()>0) nMatch++;
            if (mode) {
                cout << "jEntry = " << i << ": ";
                for (int j=0;j<matchjetlist[i].size();j++) {
                    //for (int k=0;k<matchjetlist[i][j].size();k++) 
                    cout << "{" << matchjetlist[i][j][0] << "," << matchjetlist[i][j][1] << "}";
                    cout << "\t";
                }
                cout << endl;
            }
        }
    cout << "match event: " << nMatch << "\tRatio:" << (float)nMatch/matchjetlist.size() << endl;

}
void studyHmass(string inputFile="gen2HDMbb_MZp1700_MA0300.root") {

    setNCUStyle(true);
    //get TTree from file ...
    TreeReader data(inputFile.data());
    //vector<vector<int>> genPar = genMatch_base(inputFile.data());
    vector<vector<vector<int>>> matchAk4HJetList = matchJet_4Vector(inputFile.data(),0,4);
    test(matchAk4HJetList);
    vector<vector<vector<int>>> matchAk8HJetList = matchJet_4Vector(inputFile.data(),0,8);
    test(matchAk8HJetList);
    vector<vector<vector<int>>> matchAk4A0JetList = matchJet_4Vector(inputFile.data(),1,4);
    test(matchAk4A0JetList);
    vector<vector<vector<int>>> matchAk8A0JetList = matchJet_4Vector(inputFile.data(),1,8);
    test(matchAk8A0JetList);
    TCanvas* c1 = new TCanvas("c1","c1",600,600);
    TH1F* h_hM[2][4];
    for (int j=0;j<2;j++) {
        for (int i=0;i<4;i++) {
            h_hM[j][i] = new TH1F(Form("h_AK%d_higgsM_cut%d",4+4*j,i), Form("h_AK%d_higgsM_cut%d",4+4*j,i), 50,100,150);
        }
    }
    
    for (int jEntry=0;jEntry<matchAk4HJetList.size();jEntry++) {
        data.GetEntry(jEntry);
        int nGenak4Jet = data.GetInt("ak4nGenJet");
        TClonesArray* genak4jetP4 =  (TClonesArray*) data.GetPtrTObject("ak4GenJetP4");
        for (int i=0;i<matchAk4HJetList[jEntry].size();i++) {
            TLorentzVector* thisJet = (TLorentzVector*)genak4jetP4->At(matchAk4HJetList[jEntry][i][0]);
            TLorentzVector* thatJet = (TLorentzVector*)genak4jetP4->At(matchAk4HJetList[jEntry][i][1]);
            float diJetM = (*thisJet+*thatJet).M();
            h_hM[0][0]->Fill(diJetM);
            if (diJetM>140||diJetM<110) continue;
            h_hM[0][1]->Fill(diJetM);
        }
    }
    
    // 3 jets, h1a02
    for (int jEntry=0;jEntry<matchAk8HJetList.size();jEntry++) {
        data.GetEntry(jEntry);
        float *ak8GenJetMSD = data.GetPtrFloat("ak8GenJetMSD");
        TClonesArray* genak8jetP4 =  (TClonesArray*) data.GetPtrTObject("ak8GenJetP4");
        TClonesArray* genak4jetP4 =  (TClonesArray*) data.GetPtrTObject("ak4GenJetP4");
        for (int i=0;i<matchAk8HJetList[jEntry].size();i++) {
            TLorentzVector* thisak8Jet = (TLorentzVector*)genak8jetP4->At(matchAk8HJetList[jEntry][i][0]);
            h_hM[1][0]->Fill(ak8GenJetMSD[matchAk8HJetList[jEntry][i][0]]);
            if (ak8GenJetMSD[matchAk8HJetList[jEntry][i][0]]>140 || ak8GenJetMSD[matchAk8HJetList[jEntry][i][0]]<110) continue;
            h_hM[1][1]->Fill(ak8GenJetMSD[matchAk8HJetList[jEntry][i][0]]);
            // cut A0
            for (int j=0;j<matchAk4A0JetList[jEntry].size();j++) {
                TLorentzVector* thisJet = (TLorentzVector*)genak4jetP4->At(matchAk4A0JetList[jEntry][j][0]);
                TLorentzVector* thatJet = (TLorentzVector*)genak4jetP4->At(matchAk4A0JetList[jEntry][j][1]);
                float diJetM = (*thisJet+*thatJet).M();
                if (diJetM>350||diJetM<250) continue;
                h_hM[1][2]->Fill(ak8GenJetMSD[matchAk8HJetList[jEntry][i][0]]);
                float triJetM = (*thisak8Jet+*thisJet+*thatJet).M();
                if(triJetM>1800||triJetM<1600) continue;
                h_hM[1][3]->Fill(ak8GenJetMSD[matchAk8HJetList[jEntry][i][0]]);
            }

        }
    }
    
    string pdfName = "higgsMass.pdf";
    c1->Print((pdfName+"[").data());
    for (int j=0;j<2;j++) {
        for (int i=0;i<4;i++) {
            h_hM[j][i]->Draw("hist"); 
            c1->Print(pdfName.data());
        }
    }
    c1->Print((pdfName+"]").data());

}
