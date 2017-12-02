
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
#include "../setNCUStyle.C"
#include "genMatch.C"
using namespace std;
void studyHmass(string inputFile="gen2HDMbb_MZp1700_MA0300.root") {

    setNCUStyle(true);
    //get TTree from file ...
    TreeReader data(inputFile.data());
    vector<vector<int>> genPar = genMatch_base(inputFile.data());
}
