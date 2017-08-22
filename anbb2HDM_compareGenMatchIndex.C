
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1D.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <fstream>
#include <TCanvas.h>
#include "string"
#include "./setNCUStyle.C"
#include "./anbb2HDM.h"

using namespace std;
void anbb2HDM_compareGenMatchIndex(std::string inputFile) {

    setNCUStyle(true);
    //get TTree from file ...
    TreeReader data(inputFile.data());
    
    bool upeff = false;
     
    
    int npsel = 0, npmatch = 0, ntotal = data.GetEntriesFast();
    int nHm = 0, nA0m = 0, nBothm = 0;
    for(Long64_t jEntry=1; jEntry<ntotal ;jEntry++){
        
        vector<vector<float>> selectionIdList = anbb2HDM_4bJetsselection(jEntry,inputFile);
        vector<vector<float>> genmatchIdList = anbb2HDM_4genjetMatch(jEntry,inputFile);
        if (genmatchIdList[0][0]>0) npmatch++;
        if (selectionIdList[0][0]>0) npsel++;
        if (genmatchIdList[0][0]<0||selectionIdList[0][0]<0) continue;
        bool isHm = (genmatchIdList[0][0]==selectionIdList[0][0]&&genmatchIdList[0][1]==selectionIdList[0][1])||(genmatchIdList[0][0]==selectionIdList[0][1]&&genmatchIdList[0][1]==selectionIdList[0][0]);
        bool isA0m = (genmatchIdList[0][2]==selectionIdList[0][2]&&genmatchIdList[0][3]==selectionIdList[0][3])||(genmatchIdList[0][2]==selectionIdList[0][3]&&genmatchIdList[0][3]==selectionIdList[0][2]);
        if (isHm) nHm++;
        if (isA0m) nA0m++;
        if (isHm&&isA0m) nBothm++;
    
    }
    cout << "genMatch:" << endl;
    efferr(npmatch,ntotal);
    cout << "selection:" << endl;
    efferr(npmatch,ntotal);
    cout << "Higgs match ratio: " << (float)nHm/ntotal << endl;
    cout << "A0 match ratio: " << (float)nA0m/ntotal << endl;
    cout << "Both match ratio: " << (float)nBothm/ntotal << endl;
}
