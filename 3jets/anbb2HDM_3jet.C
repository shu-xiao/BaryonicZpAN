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
void anbb2HDM_3jet(int w, std::string inputFile){
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

    TCanvas* c1 = new TCanvas("c1","",600,600);

    TH1F* h_allEvent = new TH1F("h_allEvent","h_allEvent",10,-0.5,9);
    //float bin_HT[10] = {50,100,200,300,500,700,1000,1500,2000,3000};
    //TH1F* h_HT = new TH1F("h_HT","h_HT",9,bin_HT);
    TH1F* h_HT = new TH1F("h_HT","h_HT",60,0,3000);
    
    TH1F* h_hPtAs = new TH1F("h_hPtAs","h_higgsPtAssymetry",40,0,0.8);
    TH1F* h_a0PtAs = new TH1F("h_a0PtAs","h_A0PtAssymetry",40,0,0.8);
    TH1F* h_zpPtAs = new TH1F("h_zpPtAs","h_ZpPtAssymetry",40,0,0.8);

    TH1F* h_hPtSD = new TH1F("h_hPtSD","h_higgsPtSD",30,0,0.6);
    TH1F* h_a0PtSD = new TH1F("h_a0PtSD","h_A0PtSD",30,0,0.6);
    TH1F* h_zpPtSD = new TH1F("h_zpPtSD","h_ZpPtSD",30,0,0.6);
    
    TH1F* h_hNcandi1 = new TH1F("h_hNcandi1", "h_higgs_NcandidatePairJets_1fatjet", 10,-0.5,9.5);
    TH1F* h_a0Ncandi1 = new TH1F("h_a0Ncandi1", "h_A0_NcandidatePairJets_1fatjet", 10,-0.5,9.5);
    TH1F* h_zpNcandi_h1a02 = new TH1F("h_zpNcandi_h1a02", "h_Zp_NcandidatePairJets_h1a02", 10,-0.5,9.5);
    TH1F* h_hNcandi2 = new TH1F("h_hNcandi2", "h_higgs_NcandidatePairJets_dijets", 10,-0.5,9.5);
    TH1F* h_a0Ncandi2 = new TH1F("h_a0Ncandi2", "h_A0_NcandidatePairJets_dijets", 10,-0.5,9.5);
    TH1F* h_zpNcandi_h2a01 = new TH1F("h_zpNcandi_h2a01", "h_Zp_NcandidatePairJets_h2a01", 10,-0.5,9.5);

    TH1F* h_hM = new TH1F("h_higgsM", "h_higgsM_genJet", 50,100,150);
    TH1F* h_a0M = new TH1F("h_a0M", "h_A0M_genJet", 24,A0mass.Atof()-60,A0mass.Atof()+60);
    TH1F* h_zpM = new TH1F("h_ZpM", "h_ZpM_genJet", 24,Zpmass.Atof()-120,Zpmass.Atof()+120);
    
    TH1F* h_hPt = new TH1F("h_higgsPt", "h_higgsPt_genJet", 50,0,1000);
    TH1F* h_a0Pt = new TH1F("h_a0Pt", "h_A0Pt_genJet", 40,0,800);
    
    TH1F* h_hDeltaR = new TH1F("h_HiggstobbDeltaR", "h_HiggstobbDeltaR_reco", 40,0,4);
    TH1F* h_a0DeltaR = new TH1F("h_A0tobbDeltaR", "h_A0tobbDeltaR_reco", 40,0,4);
    TH1F* h_zpDeltaR = new TH1F("h_ZptoHA0DeltaR", "h_ZptoHA0DeltaR_reco", 40,0,4);
    TH1F* h_hDeltaEta = new TH1F("h_HiggstobbDeltaEta", "h_HiggstobbDeltaEta_reco", 40,0,4);
    TH1F* h_a0DeltaEta = new TH1F("h_A0tobbDeltaEta", "h_A0tobbDeltaEta_reco", 40,0,4);
    TH1F* h_zpDeltaEta = new TH1F("h_ZptoHA0DeltaEta", "h_ZptoHA0DeltaEta_reco", 40,0,4);
    TH1F* h_hDeltaPhi = new TH1F("h_HiggstobbDeltaPhi", "h_HiggstobbDeltaPhi_reco", 32,0,3.2);
    TH1F* h_a0DeltaPhi = new TH1F("h_A0tobbDeltaPhi", "h_A0obbDeltaPhi_reco", 32,0,3.2);
    TH1F* h_zpDeltaPhi = new TH1F("h_ZptoHA0DeltaPhi", "h_ZptoHA0DeltaPhi_reco", 32,0,3.2);
    
    Int_t nPass[20]={0};
    
    //for(Long64_t jEntry=0; jEntry<100 ;jEntry++){
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        
        // broken events
        data.GetEntry(jEntry);
        if (jEntry %1000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        if (jEntry!=232) continue;
        cout << "j=" << jEntry << endl; 
        //if (jEntry>5317) cout << jEntry << endl;
        //if (jEntry==5317) continue;
        int nGenJet = (!isBG)? data.GetInt("ak4nGenJet"):data.GetInt("THINnJet");
        int nGenak8Jet = (!isBG)? data.GetInt("ak8nGenJet"):data.GetInt("FATnJet");
        float HT = data.GetFloat("HT");
        h_HT->Fill(HT);
        h_allEvent->Fill(1);
        //0. has a good vertex
        if(nGenJet<2||nGenak8Jet<1) continue;
        nPass[0]++;
        float *ak8GenJetMSD = data.GetPtrFloat("ak8GenJetMSD");
        float *ak8GenJetSDSJdR = data.GetPtrFloat("ak8GenJetSDSJdR");
        float *ak8GenJetSDSJSymm = data.GetPtrFloat("ak8GenJetSDSJSymm");
        float *ak8GenJetSDMassDrop = data.GetPtrFloat("ak8GenJetSDMassDrop");
        float *ak8GenJettau1 = data.GetPtrFloat("ak8GenJettau1");
        float *ak8GenJettau2 = data.GetPtrFloat("ak8GenJettau2");
        float *ak8GenJettau3 = data.GetPtrFloat("ak8GenJettau3");
        
        TClonesArray* genjetP4 = (!isBG)? (TClonesArray*) data.GetPtrTObject("ak4GenJetP4"):(TClonesArray*) data.GetPtrTObject("THINjetP4");
        TClonesArray* genak8jetP4 = (!isBG)? (TClonesArray*) data.GetPtrTObject("ak8GenJetP4"):(TClonesArray*) data.GetPtrTObject("THINjetP4");
        
        const int nCandidates = 200;
        vector<vector<float>> HindexList, A0indexList, ZpindexList;
        vector<vector<float>> HindexList2, A0indexList2, ZpindexList2;
        /*
        vector<float> row(5,-99);
        HindexList.assign(nCandidates,row);
        A0indexList.assign(nCandidates,row);
        ZpindexList.assign(nCandidates,row);
        HindexList2.assign(nCandidates,row);
        A0indexList2.assign(nCandidates,row);
        ZpindexList2.assign(nCandidates,row);
        */
        // reco higgs
        bool h1a02 = false, h2a01 = false;
        for(int j=0; j < nGenak8Jet; j++) {
            TLorentzVector* thisak8Jet = (TLorentzVector*)genak8jetP4->At(j);
            if(upeff && thisak8Jet->Pt()<30)continue;
            if(upeff && fabs(thisak8Jet->Eta())>2.4)continue;
            if (ak8GenJetMSD[j]<140 && ak8GenJetMSD[j]>110) {
                h1a02=true;
                HindexList2.push_back({(float)j,-1,(float)thisak8Jet->Pt()});
                //HindexList2[jList][0] = j;
                //HindexList2[jList][2] = thisak8Jet->Pt();
                //jList++;
            }
        }
        for(int ij=0, iList=0, jList=0; ij < nGenJet; ij++) {
            TLorentzVector* thisJet = (TLorentzVector*)genjetP4->At(ij);
            if(upeff && thisJet->Pt()<30)continue;
            if(upeff && fabs(thisJet->Eta())>2.4)continue;
            for (int jj=0;jj<ij;jj++) {
                TLorentzVector* thatJet = (TLorentzVector*)genjetP4->At(jj);
                if(upeff && thatJet->Pt()<30)continue;
                if(upeff && fabs(thatJet->Eta())>2.4)continue;
                float diJetM = (*thisJet+*thatJet).M();
                if (diJetM>140||diJetM<110) continue;
                HindexList.push_back({(float)ij,(float)jj,(float)(*thisJet+*thatJet).Pt()});
                /*
                HindexList[iList][0] = j;
                HindexList[iList][1] = jj;
                HindexList[iList][2] = (*thisJet+*thatJet).Pt();
                
                iList++;
                */
                h2a01=true;
                //if (abs(thisJet->Eta()-thatJet->Eta())>1.0) continue;
                //if ((*thisJet+*thatJet).Pt()>HptMax) {
                //}
            }
        } // end of outer loop jet
        cout << "HHH" << endl;
        if (HindexList.size()>0) sort(HindexList.begin(),HindexList.end(),sortListbyPt);
        if (HindexList2.size()>0) sort(HindexList2.begin(),HindexList2.end(),sortListbyPt);
        h_hNcandi1->Fill(HindexList2.size());
        h_hNcandi2->Fill(HindexList.size());
        if (HindexList.size()==0&&HindexList2.size()==0) { continue;}
        nPass[1]++;
        
        // reco A0
        if (h2a01) {
            for(int j=0; j < nGenak8Jet; j++) {                 // ak8 jet
                TLorentzVector* thisak8Jet = (TLorentzVector*)genak8jetP4->At(j);
                if(upeff && thisak8Jet->Pt()<30)continue;
                if(upeff && fabs(thisak8Jet->Eta())>2.4)continue;
                if (ak8GenJetMSD[j]<(A0mass.Atof()+50)&&ak8GenJetMSD[j]>(A0mass.Atof()-50)) {
                    //A0indexList2[jList][0] = ij;
                    //A0indexList2[jList][2] = thisak8Jet->Pt();
                    //jList++;
                    A0indexList2.push_back({(float)j,-1,(float)thisak8Jet->Pt()});
                }
            }
        }
        if (h1a02) {
            for(int ij=0, iList=0,jList=0; ij < nGenJet; ij++) {    // ak4 jet
                TLorentzVector* thisJet = (TLorentzVector*)genjetP4->At(ij);
                if(upeff && thisJet->Pt()<30)continue;
                if(upeff && fabs(thisJet->Eta())>2.4)continue;
                for (int jj=0;jj<ij;jj++) {
                    TLorentzVector* thatJet = (TLorentzVector*)genjetP4->At(jj);
                    if(upeff && thatJet->Pt()<30) continue;
                    if(upeff && fabs(thatJet->Eta())>2.4) continue;
                    float diJetM = (*thisJet+*thatJet).M();
                    if (diJetM>(A0mass.Atof()+50)||diJetM<(A0mass.Atof()-50)) continue;
                    A0indexList.push_back({(float)ij,(float)jj,(float)(*thisJet+*thisJet).Pt()});
                    
                    //A0indexList[iList][0] = j;
                    //A0indexList[iList][1] = jj;
                    //A0indexList[iList][2] = (*thisJet+*thisJet).Pt();
                    //iList++;
                    
                     
                }
            }
        } // end of outer loop jet
        cout << "a0" << endl; 
        if (A0indexList.size()>0) sort(A0indexList.begin(),A0indexList.end(),sortListbyPt);
        else h1a02=false;
        if (A0indexList2.size()>0) sort(A0indexList2.begin(),A0indexList2.end(),sortListbyPt);
        else h2a01=false;
        h_a0Ncandi1->Fill(A0indexList2.size());
        h_a0Ncandi2->Fill(A0indexList.size());
        if (A0indexList.size()==0&&A0indexList2.size()==0) {continue;}
        nPass[2]++;
        
        //reco Z'
        float zpM = -999;
        if (!(h1a02||h2a01)) continue;
        cout << "zptest" << endl;
        cout << h1a02 << h2a01 << endl;
        //if (h1a02&&h2a01) continue;
        if (h1a02) {
            for (int i=0, iList=0;i<A0indexList.size();i++) {
                for (int j=0;j<HindexList2.size();j++) {
                    TLorentzVector *diJet0 = (TLorentzVector*)genjetP4->At(A0indexList[i][0]);
                    TLorentzVector *diJet1 = (TLorentzVector*)genjetP4->At(A0indexList[i][1]);
                    TLorentzVector *fatjet = (TLorentzVector*)genak8jetP4->At(HindexList2[j][0]);
                    zpM = (*fatjet+*diJet0+*diJet1).M();
                    if (zpM>(Zpmass.Atof()+100) || zpM<(Zpmass.Atof()-100)) continue;
                    ZpindexList.push_back({(float)HindexList2[j][0],-1,(float)A0indexList[i][0],(float)A0indexList[i][1],(float)(*fatjet+*diJet0+*diJet1).Pt()});
                    
                    //ZpindexList[iList][0] = HindexList2[j][0];
                    //ZpindexList[iList][1] = -1;
                    //ZpindexList[iList][2] = A0indexList[i][0];
                    //ZpindexList[iList][3] = A0indexList[i][1];
                    //ZpindexList[iList][4] =  (*fatjet+*diJet0+*diJet1).Pt();
                    //iList++;
                    
                }
            }
        }

        cout << "ttttttttt" << endl;
        if (h2a01) {
            cout << "size: h" << HindexList.size() << " a0" << A0indexList2.size() << endl; 
            for (int i=0, iList=0;i<A0indexList2.size();i++) {
                for (int j=0;i<HindexList.size();j++) {
                    cout << "index: a0" << A0indexList2[i][0] << " h" << HindexList[j][0] << " " << HindexList[j][1] << endl;
                    TLorentzVector *fatjet = (TLorentzVector*) genak8jetP4->At(A0indexList2[i][0]);
                    TLorentzVector *diJet0 = (TLorentzVector*)genjetP4->At(HindexList[j][0]);
                    TLorentzVector *diJet1 = (TLorentzVector*)genjetP4->At(HindexList[j][1]);
                    zpM = (*fatjet+*diJet0+*diJet1).M();
                    if (zpM>(Zpmass.Atof()+100) || zpM<(Zpmass.Atof()-100)) continue;
                    ZpindexList2.push_back({(float)HindexList[j][0],(float)HindexList[j][1],(float)A0indexList2[i][0],-1,(float)(*fatjet+*diJet0+*diJet1).Pt()});
                }
            }
        }
        cout << "Zp"<< endl;
        if (ZpindexList.size()>0) sort(ZpindexList.begin(),ZpindexList.end(),sortListbyPtZp);
        else h1a02=false;
        if (ZpindexList2.size()>0) sort(ZpindexList2.begin(),ZpindexList2.end(),sortListbyPtZp);
        else h2a01=false;
        if (ZpindexList.size()==0&&ZpindexList2.size()==0) { continue;}
        nPass[3]++;
        //if (ZpindexList.size()!=1) continue;
        nPass[4]++;
        cout << "fill" << endl;
        if (h1a02) {
            h_zpNcandi_h1a02->Fill(ZpindexList.size());
            TLorentzVector* HbJet = (TLorentzVector*)genak8jetP4->At(ZpindexList[0][0]);
            h_hM->Fill(HbJet->M());
            h_hPt->Fill(HbJet->Pt());
            
            TLorentzVector* A0bJet0 = (TLorentzVector*)genjetP4->At(ZpindexList[0][2]);
            TLorentzVector* A0bJet1 = (TLorentzVector*)genjetP4->At(ZpindexList[0][3]);
            h_a0M->Fill((*A0bJet0+*A0bJet1).M());
            h_a0Pt->Fill((*A0bJet0+*A0bJet1).Pt());
            h_a0PtAs->Fill(ptAssymetry(A0bJet0,A0bJet1));
            h_a0PtSD->Fill(softDropAs(A0bJet0,A0bJet1));
            h_a0DeltaR->Fill(A0bJet0->DeltaR(*A0bJet1));
            h_a0DeltaPhi->Fill(caldePhi(A0bJet0->Phi(),A0bJet1->Phi()));
            h_a0DeltaEta->Fill(abs(A0bJet0->Eta()-A0bJet1->Eta()));
            
            TLorentzVector* A0recoJet = new TLorentzVector();
            *A0recoJet = *A0bJet0 + *A0bJet1;
            h_zpM->Fill((*HbJet+*A0bJet0+*A0bJet1).M());
            h_zpPtAs->Fill(ptAssymetry(HbJet,A0recoJet));
            h_zpPtSD->Fill(softDropAs(HbJet,A0recoJet));
            h_zpDeltaR->Fill(abs(A0recoJet->DeltaR(*HbJet)));
            h_zpDeltaPhi->Fill(caldePhi(A0recoJet->Phi(),HbJet->Phi()));
            h_zpDeltaEta->Fill(abs(HbJet->Eta()-A0recoJet->Eta()));
        }
        cout << "p1" << endl;
        cout << "list size: " << ZpindexList2.size() << endl;
        if (h2a01) {
            h_zpNcandi_h2a01->Fill(ZpindexList2.size());
            TLorentzVector* HbJet0 = (TLorentzVector*)genjetP4->At(ZpindexList2[0][0]);
            TLorentzVector* HbJet1 = (TLorentzVector*)genjetP4->At(ZpindexList2[0][1]);
            h_hM->Fill((*HbJet0+*HbJet1).M());
            h_hPt->Fill((*HbJet0+*HbJet1).Pt());
            h_hPtAs->Fill(ptAssymetry(HbJet0,HbJet1));
            h_hPtSD->Fill(softDropAs(HbJet0,HbJet1));
            h_hDeltaR->Fill(HbJet0->DeltaR(*HbJet1));
            h_hDeltaPhi->Fill(caldePhi(HbJet0->Phi(),HbJet1->Phi()));
            h_hDeltaEta->Fill(abs(HbJet0->Eta()-HbJet1->Eta()));
            
            TLorentzVector* A0bJet = (TLorentzVector*)genak8jetP4->At(ZpindexList2[0][2]);
            h_a0M->Fill(A0bJet->M());
            h_a0Pt->Fill(A0bJet->Pt());
            
            TLorentzVector* HrecoJet = new TLorentzVector();
            *HrecoJet = *HbJet0 +*HbJet1;
            h_zpM->Fill((*HbJet0+*HbJet1+*A0bJet).M());
            h_zpPtAs->Fill(ptAssymetry(HrecoJet,A0bJet));
            h_zpPtSD->Fill(softDropAs(HrecoJet,A0bJet));
            h_zpDeltaR->Fill(abs(A0bJet->DeltaR(*HrecoJet)));
            h_zpDeltaPhi->Fill(caldePhi(A0bJet->Phi(),HrecoJet->Phi()));
            h_zpDeltaEta->Fill(abs(HrecoJet->Eta()-A0bJet->Eta()));
        }
    } // end of loop over entries
    cout << "ending" << endl;    
    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    for(int i=0;i<20;i++) if(nPass[i]>0) std::cout << "nPass[" << i << "] = " << nPass[i] << std::endl;
    efferr(nPass[4],nTotal);
    return;
    if (iseffi) {
        string effifilename = Form("../njetsEffi/effi_Zpmass%s_A0mass%s_3jets.txt",Zpmass.Data(),A0mass.Data());
        fstream fp;
        fp.open(effifilename.data(), ios::out);
        fp << (float)nPass[4]/nTotal << endl;
        fp.close();
        
    }
    if (!isBG&&!iseffi) 
    {
        string pdfName = Form("anGenJet_bb2HDM_MZp%s_MA0%s.pdf",Zpmass.Data(),A0mass.Data());
        string pdfNameI = Form("anGenJet_bb2HDM_MZp%s_MA0%s.pdf(",Zpmass.Data(),A0mass.Data());
        string pdfNameF = Form("anGenJet_bb2HDM_MZp%s_MA0%s.pdf)",Zpmass.Data(),A0mass.Data());
        h_HT->Draw("hist");
        c1->Print(pdfNameI.data());
        h_zpM->Draw("hist");
        c1->Print(pdfName.data());
        h_hM->Draw("hist");
        c1->Print(pdfName.data());
        h_a0M->Draw("hist");
        c1->Print(pdfName.data());
        h_hPt->Draw("hist");
        c1->Print(pdfName.data());
        h_a0Pt->Draw("hist");
        c1->Print(pdfName.data());
        h_hPtAs->Draw("hist");
        c1->Print(pdfName.data());
        h_a0PtAs->Draw("hist");
        c1->Print(pdfName.data());
        h_zpPtAs->Draw("hist");
        c1->Print(pdfName.data());
        h_hPtSD->Draw("hist");
        c1->Print(pdfName.data());
        h_a0PtSD->Draw("hist");
        c1->Print(pdfName.data());
        h_zpPtSD->Draw("hist");
        c1->Print(pdfName.data());
        h_hDeltaR->Draw("hist");
        c1->Print(pdfName.data());
        h_a0DeltaR->Draw("hist");
        c1->Print(pdfName.data());
        h_zpDeltaR->Draw("hist");
        c1->Print(pdfName.data());
        h_hDeltaEta->Draw("hist");
        c1->Print(pdfName.data());
        h_a0DeltaEta->Draw("hist");
        c1->Print(pdfName.data());
        h_zpDeltaEta->Draw("hist");
        c1->Print(pdfName.data());
        h_hDeltaPhi->Draw("hist");
        c1->Print(pdfName.data());
        h_a0DeltaPhi->Draw("hist");
        c1->Print(pdfName.data());
        h_zpDeltaPhi->Draw("hist");
        c1->Print(pdfName.data());
        h_hNcandi1->Draw("hist");
        c1->Print(pdfName.data());
        h_hNcandi2->Draw("hist");
        c1->Print(pdfName.data());
        h_a0Ncandi1->Draw("hist");
        c1->Print(pdfName.data());
        h_a0Ncandi2->Draw("hist");
        c1->Print(pdfName.data());
        h_zpNcandi_h1a02->Draw("hist");
        c1->Print(pdfName.data());
        h_zpNcandi_h2a01->Draw("hist");
        c1->Print(pdfNameF.data());
    }
    
    string fileName;
    if (isBG) fileName = Form("QCDbg2HDMbb_%d.root",w);
    else fileName = Form("signal/bb2HDM_MZp%s_MA0%s.root",Zpmass.Data(),A0mass.Data());
    TFile* outputFile = new TFile(fileName.data(),"recreate");
    h_allEvent->Write();
    h_HT->Write();
    h_zpM->Write();
    h_hM->Write();
    h_a0M->Write();
    h_hPt->Write();
    h_a0Pt->Write();
    h_hPtAs->Write();
    h_a0PtAs->Write();
    h_zpPtAs->Write();
    h_hPtSD->Write();
    h_a0PtSD->Write();
    h_zpPtSD->Write();
    h_hDeltaR->Write();
    h_a0DeltaR->Write();
    h_zpDeltaR->Write();
    h_hDeltaEta->Write();
    h_a0DeltaEta->Write();
    h_zpDeltaEta->Write();
    h_hDeltaPhi->Write();
    h_a0DeltaPhi->Write();
    h_zpDeltaPhi->Write();
    h_hNcandi1->Write();
    h_hNcandi2->Write();
    h_a0Ncandi1->Write();
    h_a0Ncandi2->Write();
    h_zpNcandi_h1a02->Write();
    h_zpNcandi_h2a01->Write();
    outputFile->Close();

}
