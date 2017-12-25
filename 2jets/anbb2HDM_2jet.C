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
#include "../gen2HDMsample/genMatch.C"
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
void anbb2HDM_2jet(int w=0, std::string inputFile="../gen2HDMsample/gen2HDMbb_MZp1700_MA0300.root"){
    setNCUStyle(true);
    //get TTree from file ...
    TreeReader data(inputFile.data());
    vector<vector<int>> genPar = genMatch_base(inputFile.data());
    
    
    bool iseffi = false;
    bool upeff = true;
    bool isBG = false;
    
    
    TString Zpmass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*_MZp}; test=${test1%%_MA0*.root}; echo \"${test}\"",inputFile.data()));
    TString A0mass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*MA0}; test=${test1%%.root}; echo \"${test}\"",inputFile.data()));
    // set default value
    if (Zpmass.Atof()>2500||Zpmass.Atof()<500 || A0mass.Atof()>800 || A0mass.Atof()<300) {
        Zpmass = "800";
        A0mass = "300";
        isBG = true;
    }

    TCanvas* c1 = new TCanvas("c1","",600,600);

    TH1F* h_allEvent = new TH1F("h_allEvent","h_allEvent",10,-0.5,9);
    TH1F* h_HT = new TH1F("h_HT","h_HT",60,0,3000);
    
    TH1F* h_hPtAs = new TH1F("h_hPtAs","h_higgsPtAssymetry",50,0,2.5);
    TH1F* h_a0PtAs = new TH1F("h_a0PtAs","h_A0PtAssymetry",50,0,2.5);
    TH1F* h_zpPtAs = new TH1F("h_zpPtAs","h_ZpPtAssymetry",50,0,2.5);

    TH1F* h_hPtSD = new TH1F("h_hPtSD","h_higgsPtSD",30,0,0.6);
    TH1F* h_a0PtSD = new TH1F("h_a0PtSD","h_A0PtSD",30,0,0.6);
    TH1F* h_zpPtSD = new TH1F("h_zpPtSD","h_ZpPtSD",30,0,0.6);
    
    TH1F* h_hNcandi = new TH1F("h_hNcandi", "h_higgs_NcandidatePairJets", 10,-0.5,9.5);
    TH1F* h_a0Ncandi = new TH1F("h_a0Ncandi", "h_A0_NcandidatePairJets", 10,-0.5,9.5);
    TH1F* h_zpNcandi = new TH1F("h_zpNcandi", "h_Zp_NcandidatePairJets", 10,-0.5,9.5);

    TH1F* h_hM = new TH1F("h_higgsM", "h_higgsM_genJet", 50,100,150);
    TH1F* h_a0M = new TH1F("h_a0M", "h_A0M_genJet", 24,A0mass.Atof()-60,A0mass.Atof()+60);
    TH1F* h_zpM = new TH1F("h_ZpM", "h_ZpM_genJet", 24,Zpmass.Atof()-120,Zpmass.Atof()+120);
    
    TH1F* h_hPt = new TH1F("h_higgsPt", "h_higgsPt_genJet", 80,0,1600);
    TH1F* h_a0Pt = new TH1F("h_a0Pt", "h_A0Pt_genJet", 80,0,1600);
    
    TH1F* h_hDeltaR = new TH1F("h_HiggstobbDeltaR", "h_HiggstobbDeltaR_reco", 50,0,5);
    TH1F* h_a0DeltaR = new TH1F("h_A0tobbDeltaR", "h_A0tobbDeltaR_reco", 50,0,5);
    TH1F* h_zpDeltaR = new TH1F("h_ZptoHA0DeltaR", "h_ZptoHA0DeltaR_reco", 50,0,5);
    TH1F* h_hDeltaEta = new TH1F("h_HiggstobbDeltaEta", "h_HiggstobbDeltaEta_reco", 50,0,5);
    TH1F* h_a0DeltaEta = new TH1F("h_A0tobbDeltaEta", "h_A0tobbDeltaEta_reco", 40,0,4);
    TH1F* h_zpDeltaEta = new TH1F("h_ZptoHA0DeltaEta", "h_ZptoHA0DeltaEta_reco", 40,0,4);
    TH1F* h_hDeltaPhi = new TH1F("h_HiggstobbDeltaPhi", "h_HiggstobbDeltaPhi_reco", 32,0,3.2);
    TH1F* h_a0DeltaPhi = new TH1F("h_A0tobbDeltaPhi", "h_A0obbDeltaPhi_reco", 32,0,3.2);
    TH1F* h_zpDeltaPhi = new TH1F("h_ZptoHA0DeltaPhi", "h_ZptoHA0DeltaPhi_reco", 32,0,3.2);
    
    TH1F* h_hMD = new TH1F("h_hMD","h_hMD",20,0,1);
    TH1F* h_a0MD = new TH1F("h_a0MD","h_a0MD",20,0,1);
    TH1F* h_htau1 = new TH1F("h_htau1","h_htau1",20,0,1);
    TH1F* h_a0tau1 = new TH1F("h_a0tau1","h_a0tau1",20,0,1);
    TH1F* h_htau2 = new TH1F("h_htau2","h_htau2",20,0,1);
    TH1F* h_a0tau2 = new TH1F("h_a0tau2","h_a0tau2",20,0,1);
    TH1F* h_htau3 = new TH1F("h_htau3","h_htau3",20,0,1);
    TH1F* h_a0tau3 = new TH1F("h_a0tau3","h_a0tau3",20,0,1);
    TH1F* h_htau21 = new TH1F("h_htau21","h_htau21",20,0,1);
    TH1F* h_htau32 = new TH1F("h_htau32","h_htau32",20,0,1);
    TH1F* h_a0tau21 = new TH1F("h_a0tau21","h_a0tau21",20,0,1);
    TH1F* h_a0tau32 = new TH1F("h_a0tau32","h_a0tau32",20,0,1);

    
    Int_t nPass[20]={0};
    
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        
        if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        data.GetEntry(jEntry);
        float HT = data.GetFloat("HT");
        h_HT->Fill(HT);
        h_allEvent->Fill(1);
        //0. has a good vertex
        int nGenJet = (!isBG)? data.GetInt("ak8nGenJet"):data.GetInt("FATnJet");
        if(nGenJet<2) continue;
        nPass[0]++;
        
        TClonesArray* genjetP4 = (!isBG)? (TClonesArray*) data.GetPtrTObject("ak8GenJetP4"):(TClonesArray*) data.GetPtrTObject("FATjetP4");
        
        const int nCandidates = 5000;
        vector<vector<float>> HindexList, A0indexList, ZpindexList;
        vector<float> row(5,-99);
        HindexList.assign(nCandidates,row);
        A0indexList.assign(nCandidates,row);
        ZpindexList.assign(nCandidates,row);
        vector<TLorentzVector*> genHA0Par;
        TLorentzVector v1;
        if (!isBG) {
            TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
            for (int i=0;i<4;i++) genHA0Par.push_back((TLorentzVector*)genParP4->At(genPar[jEntry][i]));
        }
        else {
            for (int i=0;i<4;i++) genHA0Par.push_back(&v1);
        }
        float *ak8GenJetMSD, *ak8GenJetSDSJdR, *ak8GenJetSDSJSymm, *ak8GenJetSDMassDrop, *ak8GenJettau1, *ak8GenJettau2, *ak8GenJettau3;
        if (isBG) {
            ak8GenJetMSD = data.GetPtrFloat("FATjetSDmass");
            ak8GenJettau1 = data.GetPtrFloat("FATjetTau1");
            ak8GenJettau2 = data.GetPtrFloat("FATjetTau2");
            ak8GenJettau3 = data.GetPtrFloat("FATjetTau3");
        }
        else {
            ak8GenJetMSD = data.GetPtrFloat("ak8GenJetMSD");
            ak8GenJetSDSJdR = data.GetPtrFloat("ak8GenJetSDSJdR");
            ak8GenJetSDSJSymm = data.GetPtrFloat("ak8GenJetSDSJSymm");
            ak8GenJetSDMassDrop = data.GetPtrFloat("ak8GenJetSDMassDrop");
            ak8GenJettau1 = data.GetPtrFloat("ak8GenJettau1");
            ak8GenJettau2 = data.GetPtrFloat("ak8GenJettau2");
            ak8GenJettau3 = data.GetPtrFloat("ak8GenJettau3");
        }
        // reco higgs
        for(int ij=0, iList=0; ij < nGenJet; ij++) {
            TLorentzVector* thisJet = (TLorentzVector*)genjetP4->At(ij);
            if(upeff && thisJet->Pt()<30)continue;
            if(upeff && fabs(thisJet->Eta())>2.4)continue;
            if (isBG) {
                if(genHA0Par[0]->DeltaR(*thisJet)>0.8||genHA0Par[1]->DeltaR(*thisJet)>0.8) continue;
            }
            if (ak8GenJetMSD[ij]>140 || ak8GenJetMSD[ij]<110) continue;
            HindexList[iList][0] = ij;
            HindexList[iList][1] = ij;
            HindexList[iList][2] = thisJet->Pt();
            iList++;
            if (iList>=nCandidates) {
                cout << "Higgs index out of setting" << endl;
                return;
            }
        } // end of outer loop jet
        sort(HindexList.begin(),HindexList.end(),sortListbyPt);
        if (HindexList[0][0]<0) {h_hNcandi->Fill(0); continue;}
        nPass[1]++;
        
        for (int i=0;i<HindexList.size();i++) {
            if (HindexList[i][0]<0) {
                h_hNcandi->Fill(i);
                HindexList.erase(HindexList.begin()+i,HindexList.end());
                break;
            }
        }
        // reco A0
        for(int ij=0, iList=0; ij < nGenJet; ij++) {
            TLorentzVector* thisJet = (TLorentzVector*)genjetP4->At(ij);
            if(upeff && thisJet->Pt()<30)continue;
            if(upeff && fabs(thisJet->Eta())>2.4)continue;
            if (isBG) {
                if (genHA0Par[2]->DeltaR(*thisJet)>0.8||genHA0Par[3]->DeltaR(*thisJet)>0.8) continue;
            }
            if (ak8GenJetMSD[ij]>(A0mass.Atof()+50)||ak8GenJetMSD[ij]<(A0mass.Atof()-50)) continue;
            A0indexList[iList][0] = ij;
            A0indexList[iList][1] = ij;
            A0indexList[iList][2] = thisJet->Pt();
            iList++;
            if (iList>=nCandidates) {
                cout << "A0 index out of setting" << endl;
                return;
            }
        } // end of outer loop jet
        
        sort(A0indexList.begin(),A0indexList.end(),sortListbyPt);
        if (A0indexList[0][0]<0) {h_a0Ncandi->Fill(0); continue;}
        nPass[2]++;
        
        
        for (int i=0;i<A0indexList.size();i++) {
            if (A0indexList[i][0]<0) {
                h_a0Ncandi->Fill(i);
                A0indexList.erase(A0indexList.begin()+i,A0indexList.end());
                break;
            }
        }
        
        
        //reco Z'
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
                zpM = (*bJet0+*bJet2).M();
                if (zpM>(Zpmass.Atof()+100) || zpM<(Zpmass.Atof()-100)) continue;
                ZpindexList[iList][0] = HindexList[i][0];
                ZpindexList[iList][1] = HindexList[i][1];
                ZpindexList[iList][2] = A0indexList[j][0];
                ZpindexList[iList][3] = A0indexList[j][1];
                ZpindexList[iList][4] = (*bJet0+*bJet2).Pt();
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
        }
        if (ZpindexList.size()!=1) continue;
        nPass[4]++;
        
        TLorentzVector* HbJet0 = (TLorentzVector*)genjetP4->At(ZpindexList[0][0]);
        h_hM->Fill(ak8GenJetMSD[int(ZpindexList[0][0])]);
        h_hPt->Fill(HbJet0->Pt());
        h_htau1->Fill(ak8GenJettau1[int(ZpindexList[0][0])]);
        h_htau2->Fill(ak8GenJettau2[int(ZpindexList[0][0])]);
        h_htau3->Fill(ak8GenJettau3[int(ZpindexList[0][0])]);
        h_htau21->Fill(ak8GenJettau2[int(ZpindexList[0][0])]/ak8GenJettau1[int(ZpindexList[0][0])]);
        h_htau32->Fill(ak8GenJettau3[int(ZpindexList[0][0])]/ak8GenJettau2[int(ZpindexList[0][0])]);


        TLorentzVector* A0bJet0 = (TLorentzVector*)genjetP4->At(ZpindexList[0][2]);
        h_a0M->Fill(ak8GenJetMSD[int(ZpindexList[0][2])]);
        h_a0Pt->Fill(A0bJet0->Pt());
        h_a0tau1->Fill(ak8GenJettau1[int(ZpindexList[0][2])]);
        h_a0tau2->Fill(ak8GenJettau2[int(ZpindexList[0][2])]);
        h_a0tau3->Fill(ak8GenJettau3[int(ZpindexList[0][2])]);
        h_a0tau21->Fill(ak8GenJettau2[int(ZpindexList[0][2])]/ak8GenJettau1[int(ZpindexList[0][2])]);
        h_a0tau32->Fill(ak8GenJettau3[int(ZpindexList[0][2])]/ak8GenJettau2[int(ZpindexList[0][2])]);

        TLorentzVector* A0recoJet = new TLorentzVector(), *HrecoJet = new TLorentzVector();
        *A0recoJet = *A0bJet0;
        *HrecoJet = *HbJet0;
        h_zpM->Fill((*HbJet0+*A0bJet0).M());
        h_zpPtAs->Fill(ptAssymetry(HrecoJet,A0recoJet));
        h_zpPtSD->Fill(softDropAs(HrecoJet,A0recoJet));
        h_zpDeltaR->Fill(abs(A0recoJet->DeltaR(*HrecoJet)));
        h_zpDeltaPhi->Fill(caldePhi(A0recoJet->Phi(),HrecoJet->Phi()));
        h_zpDeltaEta->Fill(abs(HrecoJet->Eta()-A0recoJet->Eta()));
        
        if (!isBG) { 
            h_hPtAs->Fill(ak8GenJetSDSJSymm[int(ZpindexList[0][0])]);
            h_hDeltaR->Fill(ak8GenJetSDSJdR[int(ZpindexList[0][0])]);
            h_hMD->Fill(ak8GenJetSDMassDrop[int(ZpindexList[0][0])]);
            h_a0PtAs->Fill(ak8GenJetSDSJSymm[int(ZpindexList[0][2])]);
            h_a0DeltaR->Fill(ak8GenJetSDSJdR[int(ZpindexList[0][2])]);
            h_a0MD->Fill(ak8GenJetSDMassDrop[int(ZpindexList[0][2])]);
        }
    } // end of loop over entries
    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    for(int i=0;i<20;i++) if(nPass[i]>0) std::cout << "nPass[" << i << "] = " << nPass[i] << std::endl;
    efferr(nPass[4],nTotal);
    if (!isBG) 
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
        h_zpPtSD->Draw("hist");
        c1->Print(pdfName.data());
        h_hDeltaR->Draw("hist");
        c1->Print(pdfName.data());
        h_a0DeltaR->Draw("hist");
        c1->Print(pdfName.data());
        h_zpDeltaR->Draw("hist");
        c1->Print(pdfName.data());
        h_zpDeltaEta->Draw("hist");
        c1->Print(pdfName.data());
        h_zpDeltaPhi->Draw("hist");
        c1->Print(pdfName.data());
        h_hMD->Draw("hist");
        c1->Print(pdfName.data());
        h_a0MD->Draw("hist");
        c1->Print(pdfName.data());
        h_htau1->Draw("hist");
        c1->Print(pdfName.data());
        h_a0tau1->Draw("hist");
        c1->Print(pdfName.data());
        h_htau2->Draw("hist");
        c1->Print(pdfName.data());
        h_a0tau2->Draw("hist");
        c1->Print(pdfName.data());
        h_htau3->Draw("hist");
        c1->Print(pdfName.data());
        h_a0tau3->Draw("hist");
        c1->Print(pdfName.data());
        h_htau21->Draw("hist");
        c1->Print(pdfName.data());
        h_a0tau21->Draw("hist");
        c1->Print(pdfName.data());
        h_htau32->Draw("hist");
        c1->Print(pdfName.data());
        h_a0tau32->Draw("hist");
        c1->Print(pdfName.data());
        h_hNcandi->Draw("hist text0");
        c1->Print(pdfName.data());
        h_a0Ncandi->Draw("hist text0");
        c1->Print(pdfName.data());
        h_zpNcandi->Draw("hist text0");
        c1->Print(pdfNameF.data());
    }
    if (iseffi) {
        string effifilename = Form("../njetsEffi/effi_Zpmass%s_A0mass%s_2jets.txt",Zpmass.Data(),A0mass.Data());
        fstream fp;
        fp.open(effifilename.data(), ios::out);
        fp << (float)nPass[4]/nTotal << endl;
        fp.close();
        
    }
    string fileName;
    if (isBG) fileName = Form("QCDbg2HDMbb_%d.root",w);
    else fileName = Form("signal/bb2HDM_MZp%s_MA0%s.root",Zpmass.Data(),A0mass.Data());
    if (isBG) {
        TFile* outputFile = new TFile(fileName.data(),"recreate");
        h_HT->Write();
        h_zpM->Write();
        h_hM->Write();
        h_a0M->Write();
        h_zpPtAs->Write();
        h_zpPtSD->Write();
        h_hDeltaR->Write();
        h_a0DeltaR->Write();
        h_zpDeltaR->Write();
        h_zpDeltaEta->Write();
        h_zpDeltaPhi->Write();
        h_htau1->Write();
        h_a0tau1->Write();
        h_htau2->Write();
        h_a0tau2->Write();
        h_htau3->Write();
        h_a0tau3->Write();
        h_htau21->Write();
        h_a0tau21->Write();
        h_htau32->Write();
        h_a0tau32->Write();
        h_hNcandi->Write();
        h_a0Ncandi->Write();
        h_zpNcandi->Write();
    }
}
