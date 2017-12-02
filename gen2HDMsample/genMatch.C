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
        //if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        data.GetEntry(jEntry);
        
        //0. has enough jet
        bool match = false;
        int exter_H = 0, exter_A0 = 0;
        
        //TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
        int* genMomParId= data.GetPtrInt("genMomParId");
        int* genParId= data.GetPtrInt("genParId");
        int Hindex[2]={-1,-1}, A0index[2]={-1,-1};
        // match higgs
        for(int ij=0; ij < 30; ij++) {
            //TLorentzVector* thisParP4 = (TLorentzVector*)genParP4->At(ij);
            if (abs(genParId[ij])!=5) continue;
            //if(upeff && thisParP4->Pt()<30) continue;
            //if(upeff && fabs(thisParP4->Eta())>2.4) continue;
            
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
//vector <int> matchJet(TreeReader Data) {
vector <int> matchJet(string inputFile) {
    TreeReader data(inputFile.data());
    vector<vector<int>> genPar = genMatch_base(inputFile.data());
    vector<vector<int>> matchJets = {};
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        if (genPar[jEntry][4]!=1) continue;
        data.GetEntry(jEntry);
        TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
        vector<TLorentzVector*> genHA0Par;
        for (int i=0;i<4;i++) genHA0Par.push_back((TLorentzVector*)genParP4->At(genPar[jEntry][i]));
    
    }
    return {0};
}
vector<vector<vector<TLorentzVector*>>> matchJet_4Vector(string inputFile,int coneSize=4) {
    vector<vector<vector<TLorentzVector*>>> list = {};
    vector<vector<TLorentzVector*>> jet4Vect = {{}};
    TreeReader data(inputFile.data());
    vector<vector<int>> genPar = genMatch_base(inputFile.data());
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        if (genPar[jEntry][4]!=1) {
            list.push_back({{}});
            continue;
        }
        data.GetEntry(jEntry);
        vector<TLorentzVector*> genHA0Par;
        for (int i=0;i<4;i++) genHA0Par.push_back((TLorentzVector*)genParP4->At(genPar[jEntry][i]));
        TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
        int nGenJet = 0; 
        TClonesArray* genJetP4;


        if (coneSize==4) {
            nGenJet=data.GetInt("ak4nGenJet");
            genJetP4 = (TClonesArray*) data.GetPtrTObject("ak4GenJetP4");
            vector<vector<int>> matchIndex;
        
            for(int ij=0; ij < nGenJet; ij++) {
                TLorentzVector* thisJet = (TLorentzVector*)genak4jetP4->At(ij);
                if (thisJet->DeltaR(*genHA0Par[0])>0.4) continue;
                else {
                    for (int jj=0;jj<nGenJet;jj++) {
                        if (ij==jj) continue;
                        // find the repeat match
                        bool repeat = false;
                        for (int nnn=0;nnn<matchIndex.size();nnn++) {
                            if (matchIndex[nnn][0]==jj&&matchIndex[nnn][1]==ij) repeat = true
                        }
                        if (repeat) continue;
                        // search for second jet
                        TLorentzVector* thatJet = (TLorentzVector*)genak4jetP4->At(jj);
                        if (thatJet->DeltaR(*genHA0Par[1])>0.4) continue;
                        matchIndex.push_back({ij,jj});
                        jet4Vect.push_back({thisJet,thatJet});
                    }
                }
            } // end of for loop

        list.push_back(jet4Vect);
        }
        else if (coneSize==8) {
            nGenJet=data.GetInt("ak8nGenJet");
            genJetP4 = (TClonesArray*) data.GetPtrTObject("ak8GenJetP4");
        }
        else {
            cout << "coneSize = 4 or 8." << endl;
            return {{}};
        }
        
    }
}
void genMatch() {
    TreeReader data("gen2HDMbb_MZp1700_MA0300.root");
    data.GetEntry(0);
    vector<int> cc = matchJet(data);
    //TreeReader data(inputFile.data());
    //genMatch_base(data);
    //vector<vector<int>> genPar = genMatch_base("gen2HDMbb_MZp1700_MA0300.root");
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
