#ifndef ANBB2HDM
#define ANBB2HDM


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
using namespace std;
void efferr(float nsig,float ntotal,float factor=1)
{
    float eff = nsig/ntotal;
    float err = sqrt( (1-eff)*eff/ntotal);
    cout << "efficiency = " << eff*factor << " +- " << err*factor << endl;
    //ofstream myfile;
    //myfile.open();
    //myfile << eff*factor << endl;
    //myfile << err*factor;
}
float caldePhi(float phi1, float phi2) {
    float dePhi = 0;
    if (abs(phi1-phi2)>TMath::Pi()) dePhi = 2*TMath::Pi() - abs(phi1-phi2);
    else dePhi = abs(phi1-phi2);
    return dePhi;
}
bool sortListbyPt(vector<float> a, vector<float> b) {return a[2]>b[2];}
bool sortListbyPtZp(vector<float> a, vector<float> b) {return a[4]>b[4];}

vector<vector<float>> anbb2HDM_4bJetsselection(int jEntry, std::string inputFile){
//vector<vector<float>> anbb2HDM(int jEntry, TreeReader data){
    //get TTree from file ...
    TreeReader data(inputFile.data());
    
    bool upeff = false;
    //TString Zpmass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*_MZp}; test=${test1%%_MA0*.root}; echo \"${test}\"",inputFile.data()));
    //TString A0mass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*MA0}; test=${test1%%.root}; echo \"${test}\"",inputFile.data()));
     
    
    //for(Long64_t jEntry=0; jEntry<50 ;jEntry++){
    //for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        
        data.GetEntry(jEntry);
        
        //0. has enough jet
        
        TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
        int* genMomParId= data.GetPtrInt("genMomParId");
        int* genParId= data.GetPtrInt("genParId");
        vector<vector<float>> HindexList, A0indexList, ZpindexList, failList;
        vector<float> row(5,-99);
        int Hindex[2]={-1,-1}, A0index[2]={-1,-1};
        // create and initialize a 30*5 2D vector
        const int nCandidates = 30;
        HindexList.assign(nCandidates,row);
        A0indexList.assign(nCandidates,row);
        ZpindexList.assign(nCandidates,row);
        failList.assign(1,row);
        // match higgs
        for(int ij=0; ij < 30; ij++) {
            TLorentzVector* thisParP4 = (TLorentzVector*)genParP4->At(ij);
            if (abs(genParId[ij])!=5) continue;
            if(upeff && thisParP4->Pt()<30) continue;
            if(upeff && fabs(thisParP4->Eta())>2.4) continue;
            if (genMomParId[ij]==25 && Hindex[0]<0) Hindex[0] = ij;
            else if (genMomParId[ij]==25 && Hindex[1]<0) Hindex[1] = ij;
            if (genMomParId[ij]==28 && A0index[0]<0) A0index[0] = ij;
            else if (genMomParId[ij]==28 && A0index[1]<0) A0index[1] = ij;
        } // end of outer loop jet
        if (Hindex[1]<0 || Hindex[0]<0 || A0index[1]<0 || A0index[0]<0) return failList;
        
        //cout << "Event: " << jEntry << endl;
        
        TLorentzVector* HbPar0 = (TLorentzVector*)genParP4->At(Hindex[0]);
        TLorentzVector* HbPar1 = (TLorentzVector*)genParP4->At(Hindex[1]);
        TLorentzVector* A0bPar0 = (TLorentzVector*)genParP4->At(A0index[0]);
        TLorentzVector* A0bPar1 = (TLorentzVector*)genParP4->At(A0index[1]);
        //h_hM->Fill((*HbJet0+*HbJet1).M());
        
        // match Higgs
        int nGenJet = data.GetInt("ak4nGenJet");
        TClonesArray* genjetP4 = (TClonesArray*) data.GetPtrTObject("ak4GenJetP4");
        for(int ij=0, iList=0, jList=0; ij < nGenJet; ij++) {
            TLorentzVector* thisJet = (TLorentzVector*)genjetP4->At(ij);
            for (int jj=0; jj<ij;jj++) {
                TLorentzVector* thatJet = (TLorentzVector*)genjetP4->At(jj);
                //bool a0Mp=true, hMp=true;
                //float diJetM = (*thisJet+*thatJet).M();
                //if (diJetM<140&&diJetM>110) hMp=true;
                //if (diJetM<(A0mass.Atof()+50) && diJetM>(A0mass.Atof()-50)) a0Mp=true;
                if ((thisJet->DeltaR(*HbPar0)<0.4&&thatJet->DeltaR(*HbPar1)<0.4)||(thisJet->DeltaR(*HbPar1)<0.4&&thatJet->DeltaR(*HbPar0)<0.4)) {
                    HindexList[iList][0] = ij;
                    HindexList[iList][1] = jj;
                    HindexList[iList][2] = (*thisJet+*thatJet).Pt();
                    iList++;
                }
                if ((thisJet->DeltaR(*A0bPar0)<0.4&&thatJet->DeltaR(*A0bPar1)<0.4)||(thisJet->DeltaR(*A0bPar1)<0.4&&thatJet->DeltaR(*A0bPar0)<0.4)) {
                    A0indexList[jList][0] = ij;
                    A0indexList[jList][1] = jj;
                    A0indexList[jList][2] = (*thisJet+*thatJet).Pt();
                    jList++;
                }
                if (iList>=nCandidates || jList>=nCandidates) {
                    cout << "Index out of setting" << endl;
                    return failList;
                }
                
            }
        } // end of loop jet
        
        sort(HindexList.begin(),HindexList.end(),sortListbyPt);
        sort(A0indexList.begin(),A0indexList.end(),sortListbyPt);
        if (HindexList[0][0]<0)  return failList;
        // remove useless rows in 2D vector
        //cout << jEntry << endl;
        for (int i=0;i<HindexList.size();i++) {
            if (HindexList[i][0]<0) {
                HindexList.erase(HindexList.begin()+i,HindexList.end());
                break;
            }
            //cout << "Hjet" << i << " [ " << HindexList[i][0] << " , " << HindexList[i][1] << " , " << HindexList[i][2] << " ]" << endl;
        }
        
        if (A0indexList[0][0]<0) return failList;
        for (int i=0;i<A0indexList.size();i++) {
            if (A0indexList[i][0]<0) {
                A0indexList.erase(A0indexList.begin()+i,A0indexList.end());
                break;
            }
            //cout << "A0jet" << i << " [ " << A0indexList[i][0] << " , " << A0indexList[i][1] << " , " << A0indexList[i][2] << " ]" << endl;
        }
        
        //reco Z'

        //float zpM = (*HbJet0+*HbJet1+*A0bJet0+*A0bJet1).M();
        float zpM = -999;
        for (int i=0, iList=0;i<HindexList.size();i++) {
            for (int j=0;j<A0indexList.size();j++) {
                bool idCrash = false;
                for (int m=0;m<2;m++) for (int n=0;n<2;n++) if (HindexList[i][m]==A0indexList[j][n]) idCrash = true;
                if (idCrash) continue;
                TLorentzVector* bJet0 = (TLorentzVector*)genjetP4->At(HindexList[i][0]);
                TLorentzVector* bJet1 = (TLorentzVector*)genjetP4->At(HindexList[i][1]);
                TLorentzVector* bJet2 = (TLorentzVector*)genjetP4->At(A0indexList[j][0]);
                TLorentzVector* bJet3 = (TLorentzVector*)genjetP4->At(A0indexList[j][1]);
                //zpM = (*bJet0+*bJet1+*bJet2+*bJet3).M();
                //if (zpM>(Zpmass.Atof()+100) || zpM<(Zpmass.Atof()-100)) continue;
                ZpindexList[iList][0] = HindexList[i][0];
                ZpindexList[iList][1] = HindexList[i][1];
                ZpindexList[iList][2] = A0indexList[j][0];
                ZpindexList[iList][3] = A0indexList[j][1];
                ZpindexList[iList][4] = (*bJet0+*bJet1+*bJet2+*bJet3).Pt();
                iList++;
                if (iList>=nCandidates) {
                    cout << "Zp index out of setting" << endl;
                    return failList;
                }
            }
        }
        if (ZpindexList[0][0]<0) return failList;
        sort(ZpindexList.begin(),ZpindexList.end(),sortListbyPtZp);
        for (int i=0;i<ZpindexList.size();i++) {
            if (ZpindexList[i][0]<0) {
                ZpindexList.erase(ZpindexList.begin()+i,ZpindexList.end());
                break;
            }
            //cout << "jet" << i << " [ " << ZpindexList[i][0] << " , " << ZpindexList[i][1] << " , " << ZpindexList[i][2] << " ]" << endl;
        }
        
        if (ZpindexList.size()!=1) return failList;
        else return ZpindexList;

    //} // end of loop over entries


}


vector<vector<float>> anbb2HDM_4genjetMatch(int jEntry, std::string inputFile){
    //get TTree from file ...
    TreeReader data(inputFile.data());
    
    bool upeff = false;
    
    //for(Long64_t jEntry=0; jEntry<50 ;jEntry++){
    //for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        
        data.GetEntry(jEntry);
        
        //0. has enough jet
        
        TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
        int* genMomParId= data.GetPtrInt("genMomParId");
        int* genParId= data.GetPtrInt("genParId");
        vector<vector<float>> HindexList, A0indexList, ZpindexList, failList;
        vector<float> row(5,-99);
        int Hindex[2]={-1,-1}, A0index[2]={-1,-1};
        // create and initialize a 30*5 2D vector
        const int nCandidates = 30;
        HindexList.assign(nCandidates,row);
        A0indexList.assign(nCandidates,row);
        ZpindexList.assign(nCandidates,row);
        failList.assign(1,row);
        //int Hindex[2] = {-1,-1}, A0index[2] = {-1,-1};
        // match higgs
        for(int ij=0; ij < 30; ij++) {
            TLorentzVector* thisParP4 = (TLorentzVector*)genParP4->At(ij);
            if (abs(genParId[ij])!=5) continue;
            if(upeff && thisParP4->Pt()<30) continue;
            if(upeff && fabs(thisParP4->Eta())>2.4) continue;
            if (genMomParId[ij]==25 && Hindex[0]<0) Hindex[0] = ij;
            else if (genMomParId[ij]==25 && Hindex[1]<0) Hindex[1] = ij;
            if (genMomParId[ij]==28 && A0index[0]<0) A0index[0] = ij;
            else if (genMomParId[ij]==28 && A0index[1]<0) A0index[1] = ij;
        } // end of outer loop jet
        if (Hindex[1]<0 || Hindex[0]<0 || A0index[1]<0 || A0index[0]<0) return failList;
        
        //cout << "Event: " << jEntry << endl;
        
        TLorentzVector* HbPar0 = (TLorentzVector*)genParP4->At(Hindex[0]);
        TLorentzVector* HbPar1 = (TLorentzVector*)genParP4->At(Hindex[1]);
        TLorentzVector* A0bPar0 = (TLorentzVector*)genParP4->At(A0index[0]);
        TLorentzVector* A0bPar1 = (TLorentzVector*)genParP4->At(A0index[1]);
        //h_hM->Fill((*HbJet0+*HbJet1).M());
        
        // match Higgs
        int nGenJet = data.GetInt("ak4nGenJet");
        TClonesArray* genjetP4 = (TClonesArray*) data.GetPtrTObject("ak4GenJetP4");
        for(int ij=0, iList=0, jList=0; ij < nGenJet; ij++) {
            TLorentzVector* thisJet = (TLorentzVector*)genjetP4->At(ij);
            for (int jj=0; jj<ij;jj++) {
                TLorentzVector* thatJet = (TLorentzVector*)genjetP4->At(jj);
                //bool a0Mp=true, hMp=true;
                //float diJetM = (*thisJet+*thatJet).M();
                //if (diJetM<140&&diJetM>110) hMp=true;
                //if (diJetM<(A0mass.Atof()+50) && diJetM>(A0mass.Atof()-50)) a0Mp=true;
                if ((thisJet->DeltaR(*HbPar0)<0.4&&thatJet->DeltaR(*HbPar1)<0.4)||(thisJet->DeltaR(*HbPar1)<0.4&&thatJet->DeltaR(*HbPar0)<0.4)) {
                    HindexList[iList][0] = ij;
                    HindexList[iList][1] = jj;
                    HindexList[iList][2] = (*thisJet+*thatJet).Pt();
                    iList++;
                }
                if ((thisJet->DeltaR(*A0bPar0)<0.4&&thatJet->DeltaR(*A0bPar1)<0.4)||(thisJet->DeltaR(*A0bPar1)<0.4&&thatJet->DeltaR(*A0bPar0)<0.4)) {
                    A0indexList[jList][0] = ij;
                    A0indexList[jList][1] = jj;
                    A0indexList[jList][2] = (*thisJet+*thatJet).Pt();
                    jList++;
                }
                if (iList>=nCandidates || jList>=nCandidates) {
                    cout << "Index out of setting" << endl;
                    return failList;
                }
                
            }
        } // end of loop jet
        
        sort(HindexList.begin(),HindexList.end(),sortListbyPt);
        sort(A0indexList.begin(),A0indexList.end(),sortListbyPt);
        if (HindexList[0][0]<0)  return failList;
        // remove useless rows in 2D vector
        //cout << jEntry << endl;
        for (int i=0;i<HindexList.size();i++) {
            if (HindexList[i][0]<0) {
                HindexList.erase(HindexList.begin()+i,HindexList.end());
                break;
            }
            //cout << "Hjet" << i << " [ " << HindexList[i][0] << " , " << HindexList[i][1] << " , " << HindexList[i][2] << " ]" << endl;
        }
        
        if (A0indexList[0][0]<0)  return failList;
        for (int i=0;i<A0indexList.size();i++) {
            if (A0indexList[i][0]<0) {
                A0indexList.erase(A0indexList.begin()+i,A0indexList.end());
                break;
            }
            //cout << "A0jet" << i << " [ " << A0indexList[i][0] << " , " << A0indexList[i][1] << " , " << A0indexList[i][2] << " ]" << endl;
        }
        float zpM = -999;
        for (int i=0, iList=0;i<HindexList.size();i++) {
            for (int j=0;j<A0indexList.size();j++) {
                bool idCrash = false;
                for (int m=0;m<2;m++) for (int n=0;n<2;n++) if (HindexList[i][m]==A0indexList[j][n]) idCrash = true;
                if (idCrash) continue;
                TLorentzVector* bJet0 = (TLorentzVector*)genjetP4->At(HindexList[i][0]);
                TLorentzVector* bJet1 = (TLorentzVector*)genjetP4->At(HindexList[i][1]);
                TLorentzVector* bJet2 = (TLorentzVector*)genjetP4->At(A0indexList[j][0]);
                TLorentzVector* bJet3 = (TLorentzVector*)genjetP4->At(A0indexList[j][1]);
                //zpM = (*bJet0+*bJet1+*bJet2+*bJet3).M();
                //if (zpM>(Zpmass.Atof()+100) || zpM<(Zpmass.Atof()-100)) continue;
                ZpindexList[iList][0] = HindexList[i][0];
                ZpindexList[iList][1] = HindexList[i][1];
                ZpindexList[iList][2] = A0indexList[j][0];
                ZpindexList[iList][3] = A0indexList[j][1];
                ZpindexList[iList][4] = (*bJet0+*bJet1+*bJet2+*bJet3).Pt();
                iList++;
                if (iList>=nCandidates) {
                    cout << "Zp index out of setting" << endl;
                    return failList;
                }
            }
        }
        if (ZpindexList[0][0]<0)  return failList;
        sort(ZpindexList.begin(),ZpindexList.end(),sortListbyPtZp);
        for (int i=0;i<ZpindexList.size();i++) {
            if (ZpindexList[i][0]<0) {
                ZpindexList.erase(ZpindexList.begin()+i,ZpindexList.end());
                break;
            }
            //cout << "jet" << i << " [ " << ZpindexList[i][0] << " , " << ZpindexList[i][1] << " , " << ZpindexList[i][2] << " ]" << endl;
        }
        if (ZpindexList.size()!=1) return failList;
        else return ZpindexList;
        
    //} // end of loop over entries


}
#endif // ANBB2HDM
