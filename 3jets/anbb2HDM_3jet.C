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
float ptAssymetry(TLorentzVector* j1, TLorentzVector* j2) {
    float minPt = (j1->Pt()>j2->Pt())? j2->Pt():j1->Pt();
    float deR = j1->DeltaR(*j2);
    float mj = (*j1+*j2).M();
    return pow(minPt*deR/mj,2);
}
float softDropAs(TLorentzVector *j1, TLorentzVector *j2, float r0 = 0.4, float beta = 0) {
    float minPt = (j1->Pt()>j2->Pt())? j2->Pt():j1->Pt();
    return minPt/(j1->Pt()+j2->Pt())*pow(r0/j1->DeltaR(*j2),beta);
}
using namespace std;
void anbb2HDM_3jet(int w=0, std::string inputFile="../gen2HDMsample/gen2HDMbb_MZp2500_MA0300.root"){
    setNCUStyle(true);
    //get TTree from file ...
    TreeReader data(inputFile.data());
    
    
    bool upeff = true;
    bool isBG = false;
    bool iseffi = true; 
    
    TString Zpmass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*_MZp}; test=${test1%%_MA0*.root}; echo \"${test}\"",inputFile.data()));
    TString A0mass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*MA0}; test=${test1%%.root}; echo \"${test}\"",inputFile.data()));
    // set default value
    if (Zpmass.Atof()>2600||Zpmass.Atof()<600 || A0mass.Atof()>800 || A0mass.Atof()<300) {
        Zpmass = "1400";
        A0mass = "300";
        isBG = true;
    }

    TCanvas* c1 = new TCanvas("c1","c1",600,600);

    TH1F* h_allEvent = new TH1F("h_allEvent","h_allEvent",10,-0.5,9);
    //float bin_HT[10] = {50,100,200,300,500,700,1000,1500,2000,3000};
    //TH1F* h_HT = new TH1F("h_HT","h_HT",9,bin_HT);
    TH1F* h_HT = new TH1F("h_HT","h_HT",60,0,3000);
    
    TH1F* h_hPtAs_h2a01 = new TH1F("h_hPtAs_h2a01","h_higgsPtAssymetry_h2a01",50,0,2.5);
    TH1F* h_a0PtAs_h1a02 = new TH1F("h_a0PtAs_h1a02","h_A0PtAssymetry_h1a02",50,0,2.5);
    TH1F* h_zpPtAs_h1a02 = new TH1F("h_zpPtAs_h1a02","h_ZpPtAssymetry_h1a02",50,0,2.5);
    TH1F* h_zpPtAs_h2a01 = new TH1F("h_zpPtAs_h2a01","h_ZpPtAssymetry_h2a01",50,0,2.5);
    
    TH1F* h_hPtSD_h2a01 = new TH1F("h_hPtSD_h2a01","h_higgsPtSD_h2a01",30,0,0.6);
    TH1F* h_a0PtSD_h1a02 = new TH1F("h_a0PtSD_h1a02","h_A0PtSD_h1a02",30,0,0.6);
    TH1F* h_zpPtSD_h1a02 = new TH1F("h_zpPtSD_h1a02","h_ZpPtSD_h1a02",30,0,0.6);
    TH1F* h_zpPtSD_h2a01 = new TH1F("h_zpPtSD_h2a01","h_ZpPtSD_h2a01",30,0,0.6);
    
    TH1F* h_hNcandi_h1a02 = new TH1F("h_hNcandi_h1a02", "h_higgs_NcandidatePairJets_1fatjet_h1a02", 10,-0.5,9.5);
    TH1F* h_a0Ncandi_h1a02 = new TH1F("h_a0Ncandi_h1a02", "h_A0_NcandidatePairJets_dijets_h1a02", 10,-0.5,9.5);
    TH1F* h_zpNcandi_h1a02 = new TH1F("h_zpNcandi_h1a02", "h_Zp_NcandidatePairJets_h1a02", 10,-0.5,9.5);
    TH1F* h_hNcandi_h2a01 = new TH1F("h_hNcandi_h2a01", "h_higgs_NcandidatePairJets_dijets_h2a01", 10,-0.5,9.5);
    TH1F* h_a0Ncandi_h2a01 = new TH1F("h_a0Ncandi_h2a01", "h_A0_NcandidatePairJets_1fatjet_h2a01", 10,-0.5,9.5);
    TH1F* h_zpNcandi_h2a01 = new TH1F("h_zpNcandi_h2a01", "h_Zp_NcandidatePairJets_h2a01", 10,-0.5,9.5);

    TH1F* h_hM_h1a02 = new TH1F("h_higgsM_h1a02", "h_higgsM_genJet_h1a02", 50,100,150);
    TH1F* h_a0M_h1a02 = new TH1F("h_a0M_h1a02", "h_A0M_genJet_h1a02", 24,A0mass.Atof()-60,A0mass.Atof()+60);
    TH1F* h_zpM_h1a02 = new TH1F("h_ZpM_h1a02", "h_ZpM_genJet_h1a02", 24,Zpmass.Atof()-120,Zpmass.Atof()+120);
    TH1F* h_hM_h2a01 = new TH1F("h_higgsM_h2a01", "h_higgsM_genJet_h2a01", 50,100,150);
    TH1F* h_a0M_h2a01 = new TH1F("h_a0M_h2a01", "h_A0M_genJet_h2a01", 24,A0mass.Atof()-60,A0mass.Atof()+60);
    TH1F* h_zpM_h2a01 = new TH1F("h_ZpM_h2a01", "h_ZpM_genJet_h2a01", 24,Zpmass.Atof()-120,Zpmass.Atof()+120);
    
    TH1F* h_hPt_h1a02 = new TH1F("h_higgsPt_h1a02", "h_higgsPt_genJet_h1a02", 80,0,1600);
    TH1F* h_hPt_h2a01 = new TH1F("h_higgsPt_h2a01", "h_higgsPt_genJet_h2a01", 80,0,1600);
    TH1F* h_a0Pt_h1a02 = new TH1F("h_a0Pt_h1a02", "h_A0Pt_genJet_h1a02", 80,0,1600);
    TH1F* h_a0Pt_h2a01 = new TH1F("h_a0Pt_h2a01", "h_A0Pt_genJet_h2a01", 80,0,1600);
    
    TH1F* h_hDeltaR_h2a01 = new TH1F("h_HiggstobbDeltaR_h2a01", "h_HiggstobbDeltaR_reco_h2a01", 50,0,5);
    TH1F* h_a0DeltaR_h1a02 = new TH1F("h_A0tobbDeltaR_h1a02", "h_A0tobbDeltaR_reco_h1a02", 50,0,5);
    TH1F* h_zpDeltaR_h1a02 = new TH1F("h_ZptoHA0DeltaR_h1a02", "h_ZptoHA0DeltaR_reco_h1a02", 50,0,5);
    TH1F* h_zpDeltaR_h2a01 = new TH1F("h_ZptoHA0DeltaR_h2a01", "h_ZptoHA0DeltaR_reco_h2a01", 50,0,5);
    
    TH1F* h_hDeltaEta_h2a01 = new TH1F("h_HiggstobbDeltaEta_h2a01", "h_HiggstobbDeltaEta_reco_h2a01", 40,0,4);
    TH1F* h_a0DeltaEta_h1a02 = new TH1F("h_A0tobbDeltaEta_h1a02", "h_A0tobbDeltaEta_reco_h1a02", 40,0,4);
    TH1F* h_zpDeltaEta_h1a02 = new TH1F("h_ZptoHA0DeltaEta_h1a02", "h_ZptoHA0DeltaEta_reco_h1a02", 40,0,4);
    TH1F* h_zpDeltaEta_h2a01 = new TH1F("h_ZptoHA0DeltaEta_h2a01", "h_ZptoHA0DeltaEta_reco_h2a01", 40,0,4);
    
    TH1F* h_hDeltaPhi_h2a01 = new TH1F("h_HiggstobbDeltaPhi_h2a01", "h_HiggstobbDeltaPhi_reco_h2a01", 32,0,3.2);
    TH1F* h_a0DeltaPhi_h1a02 = new TH1F("h_A0tobbDeltaPhi_h1a02", "h_A0obbDeltaPhi_reco_h1a02", 32,0,3.2);
    TH1F* h_zpDeltaPhi_h1a02 = new TH1F("h_ZptoHA0DeltaPhi_h1a02", "h_ZptoHA0DeltaPhi_reco_h1a02", 32,0,3.2);
    TH1F* h_zpDeltaPhi_h2a01 = new TH1F("h_ZptoHA0DeltaPhi_h2a01", "h_ZptoHA0DeltaPhi_reco_h2a01", 32,0,3.2);
    
    
    TH1F* h_hMD = new TH1F("h_hMD","h_hMD",20,0,1);
    TH1F* h_a0MD = new TH1F("h_a0MD","h_a0MD",20,0,1);
    TH1F* h_htau1_h1a02 = new TH1F("h_htau1_h1a02","h_htau1_h1a02",20,0,1);
    TH1F* h_a0tau1_h2a01 = new TH1F("h_a0tau1_h2a01","h_a0tau1_h2a01",20,0,1);
    TH1F* h_htau2_h1a02 = new TH1F("h_htau2_h1a02","h_htau2_h1a02",20,0,1);
    TH1F* h_a0tau2_h2a01 = new TH1F("h_a0tau2_h2a01","h_a0tau2_h2a01",20,0,1);
    TH1F* h_htau3_h1a02 = new TH1F("h_htau3_h1a02","h_htau3_h1a02",20,0,1);
    TH1F* h_a0tau3_h2a01 = new TH1F("h_a0tau3_h2a01","h_a0tau3_h2a01",20,0,1);
    TH1F* h_htau21_h1a02 = new TH1F("h_htau21_h1a02","h_htau21_h1a02",20,0,1);
    TH1F* h_htau32_h1a02 = new TH1F("h_htau32_h1a02","h_htau32_h1a02",20,0,1);
    TH1F* h_a0tau21_h2a01 = new TH1F("h_a0tau21_h2a01","h_a0tau21_h2a01",20,0,1);
    TH1F* h_a0tau32_h2a01 = new TH1F("h_a0tau32_h2a01","h_a0tau32_h2a01",20,0,1);
    
    
    Int_t nPass_h1a02[20]={0};
    Int_t nPass_h2a01[20]={0};
    
    //for(Long64_t jEntry=0; jEntry<100 ;jEntry++){
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        
        // broken events
        data.GetEntry(jEntry);
        if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        int nGenak4Jet = (!isBG)? data.GetInt("ak4nGenJet"):data.GetInt("THINnJet");
        int nGenak8Jet = (!isBG)? data.GetInt("ak8nGenJet"):data.GetInt("FATnJet");
        float HT = data.GetFloat("HT");
        h_HT->Fill(HT);
        h_allEvent->Fill(1);
        //0. has a good vertex
        if(nGenak4Jet<2||nGenak8Jet<1) continue;
        nPass_h1a02[0]++;
        nPass_h2a01[0]++;
        float *ak8GenJetMSD = data.GetPtrFloat("ak8GenJetMSD");
        float *ak8GenJetSDSJdR = data.GetPtrFloat("ak8GenJetSDSJdR");
        float *ak8GenJetSDSJSymm = data.GetPtrFloat("ak8GenJetSDSJSymm");
        float *ak8GenJetSDMassDrop = data.GetPtrFloat("ak8GenJetSDMassDrop");
        float *ak8GenJettau1 = data.GetPtrFloat("ak8GenJettau1");
        float *ak8GenJettau2 = data.GetPtrFloat("ak8GenJettau2");
        float *ak8GenJettau3 = data.GetPtrFloat("ak8GenJettau3");
        
        TClonesArray* genak4jetP4 = (!isBG)? (TClonesArray*) data.GetPtrTObject("ak4GenJetP4"):(TClonesArray*) data.GetPtrTObject("THINjetP4");
        TClonesArray* genak8jetP4 = (!isBG)? (TClonesArray*) data.GetPtrTObject("ak8GenJetP4"):(TClonesArray*) data.GetPtrTObject("THINjetP4");
        
        vector<vector<float>> HindexList_h2a01, A0indexList_h1a02, ZpindexList_h1a02;
        vector<vector<float>> HindexList_h1a02, A0indexList_h2a01, ZpindexList_h2a01;
        // reco higgs
        bool h1a02 = false, h2a01 = false;
        for(int j=0; j < nGenak8Jet; j++) {
            TLorentzVector* thisak8Jet = (TLorentzVector*)genak8jetP4->At(j);
            if(upeff && thisak8Jet->Pt()<30)continue;
            if(upeff && fabs(thisak8Jet->Eta())>2.4)continue;
            if (ak8GenJetMSD[j]<140 && ak8GenJetMSD[j]>110) {
                h1a02=true;
                HindexList_h1a02.push_back({(float)j,-1,(float)thisak8Jet->Pt()});
            }
        }
        for(int ij=0; ij < nGenak4Jet; ij++) {
            TLorentzVector* thisJet = (TLorentzVector*)genak4jetP4->At(ij);
            if(upeff && thisJet->Pt()<30)continue;
            if(upeff && fabs(thisJet->Eta())>2.4)continue;
            for (int jj=0;jj<ij;jj++) {
                TLorentzVector* thatJet = (TLorentzVector*)genak4jetP4->At(jj);
                if(upeff && thatJet->Pt()<30)continue;
                if(upeff && fabs(thatJet->Eta())>2.4)continue;
                float diJetM = (*thisJet+*thatJet).M();
                if (diJetM>140||diJetM<110) continue;
                HindexList_h2a01.push_back({(float)ij,(float)jj,(float)(*thisJet+*thatJet).Pt()});
                h2a01=true;
            }
        } // end of outer loop jet
        if (HindexList_h2a01.size()>0) sort(HindexList_h2a01.begin(),HindexList_h2a01.end(),sortListbyPt);
        if (HindexList_h1a02.size()>0) sort(HindexList_h1a02.begin(),HindexList_h1a02.end(),sortListbyPt);
        h_hNcandi_h1a02->Fill(HindexList_h1a02.size());
        h_hNcandi_h2a01->Fill(HindexList_h2a01.size());
        if (h1a02) nPass_h1a02[1]++;
        if (h2a01) nPass_h2a01[1]++;
        if (HindexList_h2a01.size()==0&&HindexList_h1a02.size()==0)  continue;
        
        // reco A0
        if (h2a01) {
            for(int j=0; j < nGenak8Jet; j++) {                 // ak8 jet
                TLorentzVector* thisak8Jet = (TLorentzVector*)genak8jetP4->At(j);
                if(upeff && thisak8Jet->Pt()<30)continue;
                if(upeff && fabs(thisak8Jet->Eta())>2.4)continue;
                if (ak8GenJetMSD[j]<(A0mass.Atof()+50)&&ak8GenJetMSD[j]>(A0mass.Atof()-50)) {
                    A0indexList_h2a01.push_back({(float)j,-1,(float)thisak8Jet->Pt()});
                }
            }
        }
        if (h1a02) {
            for(int ij=0; ij < nGenak4Jet; ij++) {    // ak4 jet
                TLorentzVector* thisJet = (TLorentzVector*)genak4jetP4->At(ij);
                if(upeff && thisJet->Pt()<30)continue;
                if(upeff && fabs(thisJet->Eta())>2.4)continue;
                for (int jj=0;jj<ij;jj++) {
                    TLorentzVector* thatJet = (TLorentzVector*)genak4jetP4->At(jj);
                    if(upeff && thatJet->Pt()<30) continue;
                    if(upeff && fabs(thatJet->Eta())>2.4) continue;
                    float diJetM = (*thisJet+*thatJet).M();
                    if (diJetM>(A0mass.Atof()+50)||diJetM<(A0mass.Atof()-50)) continue;
                    A0indexList_h1a02.push_back({(float)ij,(float)jj,(float)(*thisJet+*thisJet).Pt()});
                     
                }
            }
        } // end of outer loop jet

        if (A0indexList_h1a02.size()>0) {
            sort(A0indexList_h1a02.begin(),A0indexList_h1a02.end(),sortListbyPt);
            nPass_h1a02[2]++;
        }
        else h1a02=false;
        if (A0indexList_h2a01.size()>0) {
            sort(A0indexList_h2a01.begin(),A0indexList_h2a01.end(),sortListbyPt);
            nPass_h2a01[2]++;
        }
        else h2a01=false;
        h_a0Ncandi_h2a01->Fill(A0indexList_h2a01.size());
        h_a0Ncandi_h1a02->Fill(A0indexList_h1a02.size());
        
        if (A0indexList_h1a02.size()==0&&A0indexList_h2a01.size()==0) continue;
        //reco Z'
        float zpM = -999;
        if (!(h1a02||h2a01)) continue;
        //if (h1a02&&h2a01) continue;
        if (h1a02) {
            for (int i=0;i<A0indexList_h1a02.size();i++) {
                for (int j=0;j<HindexList_h1a02.size();j++) {
                    TLorentzVector *diJet0 = (TLorentzVector*)genak4jetP4->At(A0indexList_h1a02[i][0]);
                    TLorentzVector *diJet1 = (TLorentzVector*)genak4jetP4->At(A0indexList_h1a02[i][1]);
                    TLorentzVector *fatjet = (TLorentzVector*)genak8jetP4->At(HindexList_h1a02[j][0]);
                    zpM = (*fatjet+*diJet0+*diJet1).M();
                    if (zpM>(Zpmass.Atof()+100) || zpM<(Zpmass.Atof()-100)) continue;
                    ZpindexList_h1a02.push_back({(float)HindexList_h1a02[j][0],-1,(float)A0indexList_h1a02[i][0],(float)A0indexList_h1a02[i][1],(float)(*fatjet+*diJet0+*diJet1).Pt()});
                    
                }
            }
        }

        if (h2a01) {
            for (int i=0;i<A0indexList_h2a01.size();i++) {
                for (int j=0;j<HindexList_h2a01.size();j++) {
                    TLorentzVector *fatjet = (TLorentzVector*) genak8jetP4->At(A0indexList_h2a01[i][0]);
                    TLorentzVector *diJet0 = (TLorentzVector*)genak4jetP4->At(HindexList_h2a01[j][0]);
                    TLorentzVector *diJet1 = (TLorentzVector*)genak4jetP4->At(HindexList_h2a01[j][1]);
                    zpM = (*fatjet+*diJet0+*diJet1).M();
                    if (zpM>(Zpmass.Atof()+100) || zpM<(Zpmass.Atof()-100)) continue;
                    ZpindexList_h2a01.push_back({(float)HindexList_h2a01[j][0],(float)HindexList_h2a01[j][1],(float)A0indexList_h2a01[i][0],-1,(float)(*fatjet+*diJet0+*diJet1).Pt()});
                }
            }
        }
        h_zpNcandi_h1a02->Fill(ZpindexList_h1a02.size());
        h_zpNcandi_h2a01->Fill(ZpindexList_h2a01.size());
        if (ZpindexList_h1a02.size()>0) {
            sort(ZpindexList_h1a02.begin(),ZpindexList_h1a02.end(),sortListbyPtZp);
            nPass_h1a02[3]++;
        }
        else h1a02=false;
        if (ZpindexList_h2a01.size()>0) {
            sort(ZpindexList_h2a01.begin(),ZpindexList_h2a01.end(),sortListbyPtZp);
            nPass_h2a01[3]++;
        }
        else h2a01=false;
        if (ZpindexList_h1a02.size()==0&&ZpindexList_h2a01.size()==0)  continue;
        if (ZpindexList_h1a02.size()!=1) h1a02=false;
        else nPass_h1a02[4]++;
        if (ZpindexList_h2a01.size()!=1) h2a01=false;
        else nPass_h2a01[4]++;
        // fill figures
        if (h1a02) {
            TLorentzVector* HbJet = (TLorentzVector*)genak8jetP4->At(ZpindexList_h1a02[0][0]);
            //h_hM_h1a02->Fill(HbJet->M());
            h_hM_h1a02->Fill(ak8GenJetMSD[(int)ZpindexList_h1a02[0][0]]);
            h_hPt_h1a02->Fill(HbJet->Pt());
            
            TLorentzVector* A0bJet0 = (TLorentzVector*)genak4jetP4->At(ZpindexList_h1a02[0][2]);
            TLorentzVector* A0bJet1 = (TLorentzVector*)genak4jetP4->At(ZpindexList_h1a02[0][3]);
            h_a0M_h1a02->Fill((*A0bJet0+*A0bJet1).M());
            h_a0Pt_h1a02->Fill((*A0bJet0+*A0bJet1).Pt());
            h_a0PtAs_h1a02->Fill(ptAssymetry(A0bJet0,A0bJet1));
            h_a0PtSD_h1a02->Fill(softDropAs(A0bJet0,A0bJet1));
            h_a0DeltaR_h1a02->Fill(A0bJet0->DeltaR(*A0bJet1));
            h_a0DeltaPhi_h1a02->Fill(caldePhi(A0bJet0->Phi(),A0bJet1->Phi()));
            h_a0DeltaEta_h1a02->Fill(abs(A0bJet0->Eta()-A0bJet1->Eta()));
            
            TLorentzVector* A0recoJet = new TLorentzVector();
            *A0recoJet = *A0bJet0 + *A0bJet1;
            h_zpM_h1a02->Fill((*HbJet+*A0bJet0+*A0bJet1).M());
            h_zpPtAs_h1a02->Fill(ptAssymetry(HbJet,A0recoJet));
            h_zpPtSD_h1a02->Fill(softDropAs(HbJet,A0recoJet));
            h_zpDeltaR_h1a02->Fill(abs(A0recoJet->DeltaR(*HbJet)));
            h_zpDeltaPhi_h1a02->Fill(caldePhi(A0recoJet->Phi(),HbJet->Phi()));
            h_zpDeltaEta_h1a02->Fill(abs(HbJet->Eta()-A0recoJet->Eta()));
        
            h_htau1_h1a02->Fill(ak8GenJettau1[int(ZpindexList_h1a02[0][0])]);
            h_htau2_h1a02->Fill(ak8GenJettau2[int(ZpindexList_h1a02[0][0])]);
            h_htau3_h1a02->Fill(ak8GenJettau3[int(ZpindexList_h1a02[0][0])]);
            h_htau21_h1a02->Fill(ak8GenJettau2[int(ZpindexList_h1a02[0][0])]/ak8GenJettau1[int(ZpindexList_h1a02[0][0])]);
            h_htau32_h1a02->Fill(ak8GenJettau3[int(ZpindexList_h1a02[0][0])]/ak8GenJettau2[int(ZpindexList_h1a02[0][0])]);
        }
        if (h2a01) {
            TLorentzVector* HbJet0 = (TLorentzVector*)genak4jetP4->At(ZpindexList_h2a01[0][0]);
            TLorentzVector* HbJet1 = (TLorentzVector*)genak4jetP4->At(ZpindexList_h2a01[0][1]);
            h_hM_h2a01->Fill((*HbJet0+*HbJet1).M());
            h_hPt_h2a01->Fill((*HbJet0+*HbJet1).Pt());
            h_hPtAs_h2a01->Fill(ptAssymetry(HbJet0,HbJet1));
            h_hPtSD_h2a01->Fill(softDropAs(HbJet0,HbJet1));
            h_hDeltaR_h2a01->Fill(HbJet0->DeltaR(*HbJet1));
            h_hDeltaPhi_h2a01->Fill(caldePhi(HbJet0->Phi(),HbJet1->Phi()));
            h_hDeltaEta_h2a01->Fill(abs(HbJet0->Eta()-HbJet1->Eta()));
            
            TLorentzVector* A0bJet = (TLorentzVector*)genak8jetP4->At(ZpindexList_h2a01[0][2]);
            //h_a0M_h2a01->Fill(A0bJet->M());
            h_a0M_h2a01->Fill(ak8GenJetMSD[(int)ZpindexList_h2a01[0][2]]);
            h_a0Pt_h2a01->Fill(A0bJet->Pt());
            
            TLorentzVector* HrecoJet = new TLorentzVector();
            *HrecoJet = *HbJet0 +*HbJet1;
            h_zpM_h2a01->Fill((*HbJet0+*HbJet1+*A0bJet).M());
            h_zpPtAs_h2a01->Fill(ptAssymetry(HrecoJet,A0bJet));
            h_zpPtSD_h2a01->Fill(softDropAs(HrecoJet,A0bJet));
            h_zpDeltaR_h2a01->Fill(abs(A0bJet->DeltaR(*HrecoJet)));
            h_zpDeltaPhi_h2a01->Fill(caldePhi(A0bJet->Phi(),HrecoJet->Phi()));
            h_zpDeltaEta_h2a01->Fill(abs(HrecoJet->Eta()-A0bJet->Eta()));
            
            h_a0tau1_h2a01->Fill(ak8GenJettau1[int(ZpindexList_h2a01[0][2])]);
            h_a0tau2_h2a01->Fill(ak8GenJettau2[int(ZpindexList_h2a01[0][2])]);
            h_a0tau3_h2a01->Fill(ak8GenJettau3[int(ZpindexList_h2a01[0][2])]);
            h_a0tau21_h2a01->Fill(ak8GenJettau2[int(ZpindexList_h2a01[0][2])]/ak8GenJettau1[int(ZpindexList_h2a01[0][2])]);
            h_a0tau32_h2a01->Fill(ak8GenJettau3[int(ZpindexList_h2a01[0][2])]/ak8GenJettau2[int(ZpindexList_h2a01[0][2])]);
        }
    } // end of loop over entries
    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    for(int i=0;i<20;i++) if(nPass_h1a02[i]>0) std::cout << "h1a02\tnPass[" << i << "] = " << nPass_h1a02[i] << std::endl;
    efferr(nPass_h1a02[4],nTotal);
    for(int i=0;i<20;i++) if(nPass_h1a02[i]>0) std::cout << "h2a01\tnPass[" << i << "] = " << nPass_h2a01[i] << std::endl;
    efferr(nPass_h2a01[4],nTotal);
    if (iseffi) {
        string effifilename = Form("../njetsEffi/effi_Zpmass%s_A0mass%s_3jets.txt",Zpmass.Data(),A0mass.Data());
        fstream fp;
        fp.open(effifilename.data(), ios::out);
        fp << (float)nPass_h1a02[4]/nTotal << endl;
        fp << (float)nPass_h2a01[4]/nTotal << endl;
        fp.close();
        
    }
    if (!isBG) 
    {
        string pdfName = Form("anGenJet_bb2HDM_MZp%s_MA0%s.pdf",Zpmass.Data(),A0mass.Data());
        string pdfNameI = Form("anGenJet_bb2HDM_MZp%s_MA0%s.pdf(",Zpmass.Data(),A0mass.Data());
        string pdfNameF = Form("anGenJet_bb2HDM_MZp%s_MA0%s.pdf)",Zpmass.Data(),A0mass.Data());
        h_HT->Draw("hist");
        c1->Print(pdfNameI.data());
        h_zpM_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_zpM_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_hM_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_hM_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_a0M_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_a0M_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_hPt_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_hPt_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_a0Pt_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_a0Pt_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_hPtAs_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_a0PtAs_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_zpPtAs_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_zpPtAs_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_hPtSD_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_a0PtSD_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_zpPtSD_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_zpPtSD_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_htau1_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_htau2_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_htau3_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_htau21_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_htau32_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_a0tau1_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_a0tau2_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_a0tau3_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_a0tau21_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_a0tau32_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_hDeltaR_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_a0DeltaR_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_zpDeltaR_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_zpDeltaR_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_hDeltaEta_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_a0DeltaEta_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_zpDeltaEta_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_zpDeltaEta_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_hDeltaPhi_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_a0DeltaPhi_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_zpDeltaPhi_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_zpDeltaPhi_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_hNcandi_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_hNcandi_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_a0Ncandi_h2a01->Draw("hist");
        c1->Print(pdfName.data());
        h_a0Ncandi_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_zpNcandi_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_zpNcandi_h2a01->Draw("hist");
        c1->Print(pdfNameF.data());
    }
    string fileName;
    if (isBG) fileName = Form("QCDbg2HDMbb_%d.root",w);
    else fileName = Form("signal/bb2HDM_MZp%s_MA0%s.root",Zpmass.Data(),A0mass.Data());
    if (isBG) {
        TFile* outputFile = new TFile(fileName.data(),"recreate");
        h_HT->Write();
        h_zpM_h1a02->Write();
        h_zpM_h2a01->Write();
        h_hM_h1a02->Write();
        h_hM_h2a01->Write();
        h_a0M_h1a02->Write();
        h_a0M_h2a01->Write();
        h_hPt_h1a02->Write();
        h_hPt_h2a01->Write();
        h_a0Pt_h1a02->Write();
        h_a0Pt_h2a01->Write();
        h_hPtAs_h2a01->Write();
        h_a0PtAs_h1a02->Write();
        h_zpPtAs_h1a02->Write();
        h_zpPtAs_h2a01->Write();
        h_hPtSD_h2a01->Write();
        h_a0PtSD_h1a02->Write();
        h_zpPtSD_h1a02->Write();
        h_zpPtSD_h2a01->Write();
        h_htau1_h1a02->Write();
        h_htau2_h1a02->Write();
        h_htau3_h1a02->Write();
        h_htau21_h1a02->Write();
        h_htau32_h1a02->Write();
        h_a0tau1_h2a01->Write();
        h_a0tau2_h2a01->Write();
        h_a0tau3_h2a01->Write();
        h_a0tau21_h2a01->Write();
        h_a0tau32_h2a01->Write();
        h_hDeltaR_h2a01->Write();
        h_a0DeltaR_h1a02->Write();
        h_zpDeltaR_h1a02->Write();
        h_zpDeltaR_h2a01->Write();
        h_hDeltaEta_h2a01->Write();
        h_a0DeltaEta_h1a02->Write();
        h_zpDeltaEta_h1a02->Write();
        h_zpDeltaEta_h2a01->Write();
        h_hDeltaPhi_h2a01->Write();
        h_a0DeltaPhi_h1a02->Write();
        h_zpDeltaPhi_h1a02->Write();
        h_zpDeltaPhi_h2a01->Write();
        h_hNcandi_h1a02->Write();
        h_hNcandi_h2a01->Write();
        h_a0Ncandi_h2a01->Write();
        h_a0Ncandi_h1a02->Write();
        h_zpNcandi_h1a02->Write();
        h_zpNcandi_h2a01->Write();
        outputFile->Close();
    }
}
