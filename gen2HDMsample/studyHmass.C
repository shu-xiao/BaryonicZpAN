
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
void test(vector<vector<vector<int>>> matchjetlist) {
    int nMatch = 0;
    for (int i=0;i<matchjetlist.size();i++) {
        if (matchjetlist[i].size()>0) nMatch++;
        cout << "jEntry = " << i << ": ";
        for (int j=0;j<matchjetlist[i].size();j++) {
            //for (int k=0;k<matchjetlist[i][j].size();k++) 
            cout << "{" << matchjetlist[i][j][0] << "," << matchjetlist[i][j][1] << "}";
            cout << "\t";
        }
        cout << endl;
    }
    cout << "match event: " << nMatch << "\tRatio:" << (float)nMatch/matchjetlist.size() << endl;

}
void studyHmass(string inputFile="gen2HDMbb_MZp1700_MA0300.root") {

    setNCUStyle(true);
    //get TTree from file ...
    TreeReader data(inputFile.data());
    //vector<vector<int>> genPar = genMatch_base(inputFile.data());
    vector<vector<vector<int>>> matchAk4JetList = matchJet_4Vector(inputFile.data());
    vector<vector<vector<int>>> matchAk8JetList = matchJet_4Vector(inputFile.data(),0,8);
    //test(matchAk8JetList);

}
