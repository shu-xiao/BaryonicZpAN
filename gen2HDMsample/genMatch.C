#include <iostream>
#include <stdlib.h>
#include "../untuplizer.h"
#include <vector>

using namespace std;
//vector<vector<int>> genMatch_base(TreeReader Data) {
vector<vector<int>> genMatch_base(string inputFile) {
    TreeReader data(inputFile.data());
    vector<vector<int>> genParIndexList;
    bool upeff=false;
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        data.GetEntry(jEntry);
        
        //0. has enough jet
        bool match = false;
        int exter_H = 0, exter_A0 = 0;
        
        TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
        int* genMomParId= data.GetPtrInt("genMomParId");
        int* genParId= data.GetPtrInt("genParId");
        int Hindex[2]={-1,-1}, A0index[2]={-1,-1};
        // match higgs
        for(int ij=0; ij < 30; ij++) {
            TLorentzVector* thisParP4 = (TLorentzVector*)genParP4->At(ij);
            if (abs(genParId[ij])!=5) continue;
            if(upeff && thisParP4->Pt()<30) continue;
            if(upeff && fabs(thisParP4->Eta())>2.4) continue;
            
            if (genMomParId[ij]==25 && Hindex[0]<0) Hindex[0] = ij;
            else if (genMomParId[ij]==25 && Hindex[1]<0) Hindex[1] = ij;
            else if (genMomParId[ij]==25) exter_H++;
            
            if (genMomParId[ij]==28 && A0index[0]<0) A0index[0] = ij;
            else if (genMomParId[ij]==28 && A0index[1]<0) A0index[1] = ij;
            else if (genMomParId[ij]==28) exter_A0++;
        } // end of outer loop jet
        if (Hindex[1]>=0&&A0index[1]>=0) match=true;
        genParIndexList.push_back({Hindex[0],Hindex[1],A0index[0],A0index[1],match,exter_H,exter_A0});
    }
    return genParIndexList;
}
void genMatch() {
    TreeReader data("gen2HDMbb_MZp1700_MA0300.root");
    //TreeReader data(inputFile.data());
    //genMatch_base(data);
    vector<vector<int>> genPar = genMatch_base("gen2HDMbb_MZp1700_MA0300.root");
    /*
    cout << a.size() << endl;
    for (int i=0;i<10000;i++) {
        for (int j=0;j<a[i].size();j++) {
            cout << a[i][j];
            if (j<4) cout << "\t";
            else cout << " ";
        }
        cout << endl;
    }
    */
}
