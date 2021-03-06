
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
#include "setNCUStyle.C"
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

using namespace std;
void anbb2HDM_4genjetMatch(int w, std::string inputFile){
    setNCUStyle(true);
    //get TTree from file ...
    TreeReader data(inputFile.data());
    
    bool upeff = false;
    TString Zpmass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*_MZp}; test=${test1%%_MA0*.root}; echo \"${test}\"",inputFile.data()));
    TString A0mass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*MA0}; test=${test1%%.root}; echo \"${test}\"",inputFile.data()));

    TCanvas* c1 = new TCanvas("c1","",600,600);

    TH1F* h_hpt = new TH1F("h_higgsPt", "h_higgsPt", 50,0,1000);
    TH1F* h_a0pt = new TH1F("h_a0Pt", "h_a0Pt", 50,0,1000);
    
    TH1F* h_hNcandi = new TH1F("h_hNcandi", "h_higgs_NcandidatePairJets", 10,-0.5,9.5);
    TH1F* h_a0Ncandi = new TH1F("h_a0Ncandi", "h_A0_NcandidatePairJets", 10,-0.5,9.5);
    TH1F* h_zpNcandi = new TH1F("h_zpNcandi", "h_Zp_NcandidatePairJets", 10,-0.5,9.5);

    TH1F* h_hM = new TH1F("h_higgsM", "h_higgsM_genJet", 50,70,180);
    TH1F* h_a0M = new TH1F("h_a0M", "h_A0M_genJet", 40,A0mass.Atof()-100,A0mass.Atof()+100);
    TH1F* h_zpM = new TH1F("h_ZpM", "h_ZpM_genJet", 60,Zpmass.Atof()-300,Zpmass.Atof()+300);
    
    TH1F* h_hPtAs = new TH1F("h_HiggsPtAs","h_HiggsPtAsymmetry",50,0,0.5);
    TH1F* h_A0PtAs = new TH1F("h_A0PtAs","h_A0PtAsymmetry",50,0,0.5);
    TH1F* h_hPtSDAs = new TH1F("h_HiggsSDAs","h_HiggsSoftDropConditionAsymmetry",50,0,0.5);
    TH1F* h_a0PtSDAs = new TH1F("h_A0SDAs","h_A0SoftDropConditionAsymmetry",50,0,0.5);
    TH1F* h_hGamma = new TH1F("h_HiggsJetGamma", "h_HiggsJetGamma", 50,0,10);
    TH1F* h_a0Gamma = new TH1F("h_A0JetGamma", "h_A0JetGamma", 50,0,10);
    
    TH1F* h_hDeltaR = new TH1F("h_HiggstobbdeltaR", "h_HiggstobbDeltaR_genmatch", 50,0,5);
    TH1F* h_hDeltaPhi = new TH1F("h_HiggstobbdeltaPhi", "h_HiggstobbDeltaPhi_genmatch", 32,0,3.2);
    TH1F* h_hDeltaEta = new TH1F("h_HiggstobbdeltaEta", "h_HiggstobbDeltaEta_genmatch", 40,0,4);
    TH1F* h_a0DeltaR = new TH1F("h_A0tobbdeltaR", "h_A0tobbDeltaR_genmatch", 50,0,5);
    TH1F* h_a0DeltaPhi = new TH1F("h_A0tobbdeltaPhi", "h_A0tobbDeltaPhi_genmatch", 32,0,3.2);
    TH1F* h_a0DeltaEta = new TH1F("h_A0tobbdeltaEta", "h_A0tobbDeltaEta_genmatch", 40,0,4);
    TH1F* h_zpDeltaR = new TH1F("h_ZptoA0HdeltaR", "h_ZptoA0HDeltaR_genmatch", 50,0,5);
     
    TH1F* h_2candiHdePt = new TH1F("h_2candiHdePt", "h_2candiHdePt", 50,0,300);
    TH1F* h_2candiHdeDeltaR = new TH1F("h_2candiHdeDeltaR", "h_2candiHdeDeltaR", 40,0,0.4);
    TH1F* h_2candiA0dePt = new TH1F("h_2candiA0dePt", "h_2candiA0dePt", 50,0,300);
    TH1F* h_2candiA0deDeltaR = new TH1F("h_2candiA0deDeltaR", "h_2candiA0deDeltaR", 40,0,0.4);
    TH1F* h_2candiHmindePt = new TH1F("h_2candiHmindePt", "h_2candiHmindePt", 50,0,300);
    TH1F* h_2candiHmindeDeltaR = new TH1F("h_2candiHmindeDeltaR", "h_2candiHmindeDeltaR", 40,0,0.4);
    TH1F* h_2candiA0mindePt = new TH1F("h_2candiA0mindePt", "h_2candiA0mindePt", 50,0,300);
    TH1F* h_2candiA0mindeDeltaR = new TH1F("h_2candiA0mindeDeltaR", "h_2candiA0mindeDeltaR", 40,0,0.4);
    Int_t nPass[20]={0};
    
    //for(Long64_t jEntry=0; jEntry<50 ;jEntry++){
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        
        if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        data.GetEntry(jEntry);
        
        //0. has enough jet
        
        TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
        int* genMomParId= data.GetPtrInt("genMomParId");
        int* genParId= data.GetPtrInt("genParId");
        vector<vector<float>> HindexList, A0indexList, ZpindexList, ak8ZpindexList;
        vector<vector<float>> ak8HindexList, ak8A0indexList;
        vector<float> row(5,-99);
        int Hindex[2]={-1,-1}, A0index[2]={-1,-1};
        // create and initialize a 30*5 2D vector
        const int nCandidates = 30;
        HindexList.assign(nCandidates,row);
        A0indexList.assign(nCandidates,row);
        ZpindexList.assign(nCandidates,row);
        ak8HindexList.assign(nCandidates,row);
        ak8A0indexList.assign(nCandidates,row);
        ak8ZpindexList.assign(nCandidates,row);
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
        if (Hindex[1]<0 || Hindex[0]<0 || A0index[1]<0 || A0index[0]<0) continue;
        nPass[0]++;
        
        //cout << "Event: " << jEntry << endl;
        
        TLorentzVector* HbPar0 = (TLorentzVector*)genParP4->At(Hindex[0]);
        TLorentzVector* HbPar1 = (TLorentzVector*)genParP4->At(Hindex[1]);
        TLorentzVector* A0bPar0 = (TLorentzVector*)genParP4->At(A0index[0]);
        TLorentzVector* A0bPar1 = (TLorentzVector*)genParP4->At(A0index[1]);
        //h_hM->Fill((*HbJet0+*HbJet1).M());
        
        // match Higgs
        int nGenJet = data.GetInt("ak4nGenJet");
        TClonesArray* genjetP4 = (TClonesArray*) data.GetPtrTObject("ak4GenJetP4");
        for(int ij=0, iList=0, jList=0, il=0,jl=0; ij < nGenJet; ij++) {
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
                    return;
                }
                
            }
        } // end of loop jet

        int nGenak8Jet = data.GetInt("ak8nGenJet");
        TClonesArray* genak8jetP4 = (TClonesArray*) data.GetPtrTObject("ak8GenJetP4");
        for (int il=0,jl=0, ii=0; ii<nGenak8Jet;ii++) {
            TLorentzVector* thisJet = (TLorentzVector*)genak8jetP4->At(ij);
            if ((thisJet->DeltaR(*A0bPar0)<0.8&&thisJet->DeltaR(*A0bPar1)<0.8))
            {
                ak8A0indexList[il][0] = ii;
                ak8A0indexList[il][2] = thisJet->Pt();
                il++;
            }
            if((thisJet->DeltaR(*HbPar0)<0.8&&thitJet->DeltaR(*HbPar1)<0.8)) {
                ak8HindexList[jl][0] = ii;
                ak8HindexList[jl][2] = thisJet->Pt();
                jl++;
            }
            if (il>=nCandidates || jl>=nCandidates) {
                cout << "Index out of setting" << endl;
                return;
            }
        
        }
        sort(ak8HindexList.begin(),ak8HindexList.end(),sortListbyPt);
        sort(ak8A0indexList.begin(),ak8A0indexList.end(),sortListbyPt);
        sort(HindexList.begin(),HindexList.end(),sortListbyPt);
        sort(A0indexList.begin(),A0indexList.end(),sortListbyPt);
        if (HindexList[0][0]<0) {h_hNcandi->Fill(0); continue;}
        nPass[1]++;
        // remove useless rows in 2D vector
        //cout << jEntry << endl;
        for (int i=0;i<HindexList.size();i++) {
            if (HindexList[i][0]<0) {
                h_hNcandi->Fill(i);
                HindexList.erase(HindexList.begin()+i,HindexList.end());
                break;
            }
            //cout << "Hjet" << i << " [ " << HindexList[i][0] << " , " << HindexList[i][1] << " , " << HindexList[i][2] << " ]" << endl;
        }
        
        if (A0indexList[0][0]<0) {h_a0Ncandi->Fill(0); continue;}
        nPass[2]++;
        for (int i=0;i<A0indexList.size();i++) {
            if (A0indexList[i][0]<0) {
                h_a0Ncandi->Fill(i);
                A0indexList.erase(A0indexList.begin()+i,A0indexList.end());
                break;
            }
            //cout << "A0jet" << i << " [ " << A0indexList[i][0] << " , " << A0indexList[i][1] << " , " << A0indexList[i][2] << " ]" << endl;
        }
        
        //TLorentzVector* A0bJet0 = (TLorentzVector*)genjetP4->At(A0indexList[0][0]);
        //TLorentzVector* A0bJet1 = (TLorentzVector*)genjetP4->At(A0indexList[0][1]);
        //h_a0M->Fill((*A0bJet0+*A0bJet1).M());
        
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
                    return;
                }
            }
        }
        if (ZpindexList[0][0]<0) {h_zpNcandi->Fill(0); continue;}
        nPass[3]++;
        sort(ZpindexList.begin(),ZpindexList.end(),sortListbyPtZp);
        for (int i=0;i<ZpindexList.size();i++) {
            if (ZpindexList[i][0]<0) {
                h_zpNcandi->Fill(i);
                ZpindexList.erase(ZpindexList.begin()+i,ZpindexList.end());
                break;
            }
            //cout << "jet" << i << " [ " << ZpindexList[i][0] << " , " << ZpindexList[i][1] << " , " << ZpindexList[i][2] << " ]" << endl;
        }
        /*
        if (ZpindexList.size()!=2) continue;   
        vector <int> n;
        for (int i=0;i<4;i++) {
            if (ZpindexList[0][i]!=ZpindexList[1][i]) n.push_back(i);
        }
        if (n.size()<2) continue;
        TLorentzVector* j1 = (TLorentzVector*)genjetP4->At(ZpindexList[1][n[0]]);
        TLorentzVector* j2 = (TLorentzVector*)genjetP4->At(ZpindexList[1][n[1]]);
        TLorentzVector* jj1 = (TLorentzVector*)genjetP4->At(ZpindexList[0][n[0]]);
        TLorentzVector* jj2 = (TLorentzVector*)genjetP4->At(ZpindexList[0][n[1]]);
        float R1 = (j1->DeltaR(*jj1)>j1->DeltaR(*jj2))? j1->DeltaR(*jj2):j1->DeltaR(*jj1);
        float R2 = (j2->DeltaR(*jj2)>j2->DeltaR(*jj1))? j2->DeltaR(*jj1):j2->DeltaR(*jj2);
        cout << "Par dR: " << HbPar0->DeltaR(*HbPar1) << endl;;
        cout << "beC:" << j1->DeltaR(*jj1) << "\t" << j2->DeltaR(*jj2) << endl;
        cout << "R1: " << R1 << "\tR2: " << R2 << endl;
        cout << endl;
        
        for (int i=0;i<4;i++) {
            if (ZpindexList[0][i]!=ZpindexList[1][i]&&ZpindexList[0][i+1]!=ZpindexList[1][i-1]) {
                TLorentzVector* j1 = (TLorentzVector*)genjetP4->At(ZpindexList[1][i]);
                TLorentzVector* j2 = (TLorentzVector*)genjetP4->At(ZpindexList[0][i]);
                //float R1 = 
                cout << "jEntry: " << jEntry << "\tindex: " << i << "\tdeltaR: " << j1->DeltaR(*j2) <<endl;
                //j1->Print();
                //j2->Print();
                
            }
        }
        */


        TLorentzVector* HbJet0 = (TLorentzVector*)genjetP4->At(ZpindexList[0][0]);
        TLorentzVector* HbJet1 = (TLorentzVector*)genjetP4->At(ZpindexList[0][1]);
        TLorentzVector* HJet = new TLorentzVector();
        *HJet = *HbJet0 + *HbJet1;
        h_hpt->Fill(HJet->Pt());
        h_hDeltaR->Fill(HbJet0->DeltaR(*HbJet1));
        h_hDeltaEta->Fill(abs(HbJet0->Eta()-HbJet1->Eta()));
        h_hDeltaPhi->Fill(caldePhi(HbJet0->Phi(),HbJet1->Phi()));
        h_hM->Fill(HJet->M());
        h_hGamma->Fill(HJet->Gamma());
        float minHbPt = (HbJet0->Pt()>HbJet1->Pt()) ? HbJet1->Pt() : HbJet0->Pt();
        h_hPtAs->Fill(pow(minHbPt,2)/pow(HJet->M(),2)*pow(HbJet0->DeltaR(*HbJet1),2));
        
        TLorentzVector* A0bJet0 = (TLorentzVector*)genjetP4->At(ZpindexList[0][2]);
        TLorentzVector* A0bJet1 = (TLorentzVector*)genjetP4->At(ZpindexList[0][3]);
        TLorentzVector* A0Jet = new TLorentzVector();
        *A0Jet = *A0bJet0 + *A0bJet1;
        h_a0pt->Fill(A0Jet->Pt());
        h_a0DeltaR->Fill(A0bJet0->DeltaR(*A0bJet1));
        h_a0DeltaEta->Fill(abs(A0bJet0->Eta()-A0bJet1->Eta()));
        h_a0DeltaPhi->Fill(caldePhi(A0bJet0->Phi(),A0bJet1->Phi()));
        h_a0M->Fill(A0Jet->M());
        h_a0Gamma->Fill(A0Jet->Gamma());
        
        h_zpM->Fill((*HbJet0+*HbJet1+*A0bJet0+*A0bJet1).M());
        h_zpDeltaR->Fill(A0Jet->DeltaR(*HJet));
   
        TLorentzVector *bJet[2][4];
        if (ZpindexList.size()==2) {
            float mindeR[4]={9,9,9,9} , mindePt[4] = {999,999,999,999};
            for (int i=0;i<2;i++) {
                for (int j=0;j<4;j++) bJet[i][j] = (TLorentzVector*)genjetP4->At(ZpindexList[i][j]);
            }
            for (int i=0;i<2;i++) {
                h_2candiHdePt->Fill(abs(bJet[i][0]->Pt()-HbPar0->Pt()));
                h_2candiHdePt->Fill(abs(bJet[i][1]->Pt()-HbPar1->Pt()));
                h_2candiA0dePt->Fill(abs(bJet[i][2]->Pt()-A0bPar0->Pt()));
                h_2candiA0dePt->Fill(abs(bJet[i][3]->Pt()-A0bPar1->Pt()));
                h_2candiHdeDeltaR->Fill(bJet[i][0]->DeltaR(*HbPar0));
                h_2candiHdeDeltaR->Fill(bJet[i][1]->DeltaR(*HbPar1));
                h_2candiA0deDeltaR->Fill(bJet[i][2]->DeltaR(*A0bPar0));
                h_2candiA0deDeltaR->Fill(bJet[i][3]->DeltaR(*A0bPar1));
            }
            mindeR[0] = (bJet[0][0]->DeltaR(*HbPar0) > bJet[1][0]->DeltaR(*HbPar0))? bJet[1][0]->DeltaR(*HbPar0): bJet[0][0]->DeltaR(*HbPar0);
            mindeR[1] = (bJet[0][1]->DeltaR(*HbPar1) > bJet[1][1]->DeltaR(*HbPar1))? bJet[1][1]->DeltaR(*HbPar1): bJet[0][1]->DeltaR(*HbPar1);
            mindeR[2] = (bJet[0][2]->DeltaR(*A0bPar0) > bJet[1][2]->DeltaR(*A0bPar0))? bJet[1][2]->DeltaR(*A0bPar0): bJet[0][2]->DeltaR(*A0bPar0);
            mindeR[3] = (bJet[0][3]->DeltaR(*A0bPar1) > bJet[1][3]->DeltaR(*A0bPar1))? bJet[1][3]->DeltaR(*A0bPar1): bJet[0][3]->DeltaR(*A0bPar1);
            
            mindePt[0] = (abs(bJet[0][0]->Pt()-HbPar0->Pt()) > abs(bJet[1][0]->Pt()-HbPar0->Pt())) ? abs(bJet[1][0]->Pt()-HbPar0->Pt()): abs(bJet[0][0]->Pt()-HbPar0->Pt());
            mindePt[1] = (abs(bJet[0][1]->Pt()-HbPar1->Pt()) > abs(bJet[1][1]->Pt()-HbPar1->Pt())) ? abs(bJet[1][1]->Pt()-HbPar1->Pt()): abs(bJet[0][1]->Pt()-HbPar0->Pt());
            mindePt[2] = (abs(bJet[0][2]->Pt()-A0bPar0->Pt()) > abs(bJet[1][2]->Pt()-A0bPar0->Pt())) ? abs(bJet[1][2]->Pt()-A0bPar0->Pt()): abs(bJet[0][2]->Pt()-A0bPar0->Pt());
            mindePt[3] = (abs(bJet[0][3]->Pt()-A0bPar1->Pt()) > abs(bJet[1][3]->Pt()-A0bPar1->Pt())) ? abs(bJet[1][3]->Pt()-A0bPar1->Pt()): abs(bJet[0][3]->Pt()-A0bPar1->Pt());
            h_2candiHmindeDeltaR->Fill(mindeR[0]);
            h_2candiHmindeDeltaR->Fill(mindeR[1]);
            h_2candiHmindePt->Fill(mindePt[0]);
            h_2candiHmindePt->Fill(mindePt[1]);
            h_2candiA0mindeDeltaR->Fill(mindeR[2]);
            h_2candiA0mindeDeltaR->Fill(mindeR[3]);
            h_2candiA0mindePt->Fill(mindePt[2]);
            h_2candiA0mindePt->Fill(mindePt[3]);
        }
        bool isH=false, isA0=false,isZp=false;
        for (int i=0;i<ZpindexList.size();i++) {
            TLorentzVector* bJet0 = (TLorentzVector*)genjetP4->At(ZpindexList[i][0]);
            TLorentzVector* bJet1 = (TLorentzVector*)genjetP4->At(ZpindexList[i][1]);
            TLorentzVector* bJet2 = (TLorentzVector*)genjetP4->At(ZpindexList[i][2]);
            TLorentzVector* bJet3 = (TLorentzVector*)genjetP4->At(ZpindexList[i][3]);
            float Hmass = (*bJet0+*bJet1).M();
            float A0mass = (*bJet2+*bJet3).M();
            float Zpmass = (*bJet0+*bJet1+*bJet2+*bJet3).M();
            if (Hmass>110&&Hmass<140) isH=true;
            if (A0mass<350&&A0mass>250) isA0=true;
            if (Zpmass<900&&Zpmass>700) isZp=true;
        }
        if (!isH) continue;
        nPass[4]++;
        if (!isA0) continue;
        nPass[5]++;
        if (!isZp) continue;
        nPass[6]++;
        if (ZpindexList.size()!=1) continue;
        nPass[7]++;
    } // end of loop over entries
    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    for(int i=0;i<20;i++) if(nPass[i]>0) std::cout << "nPass[" << i << "] = " << nPass[i] << std::endl;
    efferr(nPass[7],nTotal);
    string pdfName = Form("genMatch_bb2HDM_MZp%s_MA0%s.pdf",Zpmass.Data(),A0mass.Data());
    string pdfNameI = Form("genMatch_bb2HDM_MZp%s_MA0%s.pdf(",Zpmass.Data(),A0mass.Data());
    string pdfNameF = Form("genMatch_bb2HDM_MZp%s_MA0%s.pdf)",Zpmass.Data(),A0mass.Data());
    h_zpM->Draw("hist");
    c1->Print(pdfNameI.data());
    h_hM->Draw("hist");
    c1->Print(pdfName.data());
    h_a0M->Draw("hist");
    c1->Print(pdfName.data());
    h_hpt->Draw("hist");
    c1->Print(pdfName.data());
    h_a0pt->Draw("hist");
    c1->Print(pdfName.data());
    h_hDeltaR->Draw("hist");
    c1->Print(pdfName.data());
    h_a0DeltaR->Draw("hist");
    c1->Print(pdfName.data());
    h_hDeltaEta->Draw("hist");
    c1->Print(pdfName.data());
    h_a0DeltaEta->Draw("hist");
    c1->Print(pdfName.data());
    h_hDeltaPhi->Draw("hist");
    c1->Print(pdfName.data());
    h_a0DeltaPhi->Draw("hist");
    c1->Print(pdfName.data());
    h_hGamma->Draw("hist");
    c1->Print(pdfName.data());
    h_a0Gamma->Draw("hist");
    c1->Print(pdfName.data());
    h_hNcandi->Draw("hist");
    c1->Print(pdfName.data());
    h_a0Ncandi->Draw("hist");
    c1->Print(pdfName.data());
    h_zpNcandi->Draw("hist");
    c1->Print(pdfName.data());
    h_2candiHdePt->Draw("hist");
    c1->Print(pdfName.data());
    h_2candiHmindePt->Draw("hist");
    c1->Print(pdfName.data());
    h_2candiA0dePt->Draw("hist");
    c1->Print(pdfName.data());
    h_2candiA0mindePt->Draw("hist");
    c1->Print(pdfName.data());
    h_2candiHdeDeltaR->Draw("hist");
    c1->Print(pdfName.data());
    h_2candiHmindeDeltaR->Draw("hist");
    c1->Print(pdfName.data());
    h_2candiA0deDeltaR->Draw("hist");
    c1->Print(pdfName.data());
    h_2candiA0mindeDeltaR->Draw("hist");
    c1->Print(pdfNameF.data());
    /* 
    string fileName = Form("bb2HDM_%d.root",w);
    //string fileName = Form("QCDbg2HDMbb_%d.root",w);
    TFile* outputFile = new TFile(fileName.data(),"recreate");
    h_hpt->Write();
    h_hM->Write();
    h_hDeltaR->Write();
    h_hDeltaPhi->Write();
    h_hDeltaEta->Write();
    h_hpassCISt0->Write();
    h_hpassCISm0->Write();
    h_hpassCISl0->Write();
    h_hpassCISt1->Write();
    h_hpassCISm1->Write();
    h_hpassCISl1->Write();
    h_passCISt->Write();
    h_passCISm->Write();
    h_passCISl->Write();
    outputFile->Close();
*/


}
