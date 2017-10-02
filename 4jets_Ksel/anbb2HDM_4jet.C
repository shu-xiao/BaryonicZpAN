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
void anbb2HDM_4jet(int w, std::string inputFile){
    setNCUStyle(true);
    //get TTree from file ...
    TreeReader data(inputFile.data());
    
    
    bool upeff = true;
    bool isBG = false;
    
    
    TString Zpmass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*_MZp}; test=${test1%%_MA0*.root}; echo \"${test}\"",inputFile.data()));
    TString A0mass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*MA0}; test=${test1%%.root}; echo \"${test}\"",inputFile.data()));
    // set default value
    if (Zpmass.Atof()>2500||Zpmass.Atof()<800 || A0mass.Atof()>800 || A0mass.Atof()<300) {
        Zpmass = "800";
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
    
    TH1F* h_hNcandi = new TH1F("h_hNcandi", "h_higgs_NcandidatePairJets", 10,-0.5,9.5);
    TH1F* h_a0Ncandi = new TH1F("h_a0Ncandi", "h_A0_NcandidatePairJets", 10,-0.5,9.5);
    TH1F* h_zpNcandi = new TH1F("h_zpNcandi", "h_Zp_NcandidatePairJets", 10,-0.5,9.5);

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
        if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        data.GetEntry(jEntry);
        float HT = data.GetFloat("HT");
        h_HT->Fill(HT);
        h_allEvent->Fill(1);
        //0. has a good vertex
        int nGenJet = (!isBG)? data.GetInt("ak4nGenJet"):data.GetInt("THINnJet");
        if(nGenJet<4) continue;
        nPass[0]++;
        
        TClonesArray* genjetP4 = (!isBG)? (TClonesArray*) data.GetPtrTObject("ak4GenJetP4"):(TClonesArray*) data.GetPtrTObject("THINjetP4");
        
        const int nCandidates = 30;
        vector<vector<float>> HindexList, A0indexList, ZpindexList;
        vector<float> row(5,-99);
        HindexList.assign(nCandidates,row);
        A0indexList.assign(nCandidates,row);
        ZpindexList.assign(nCandidates,row);
        // reco higgs
        for(int ij=0, iList=0; ij < nGenJet; ij++) {
            TLorentzVector* thisJet = (TLorentzVector*)genjetP4->At(ij);
            if(upeff && thisJet->Pt()<30)continue;
            if(upeff && fabs(thisJet->Eta())>2.4)continue;
            for (int jj=0;jj<ij;jj++) {
                TLorentzVector* thatJet = (TLorentzVector*)genjetP4->At(jj);
                if(upeff && thatJet->Pt()<30)continue;
                if(upeff && fabs(thatJet->Eta())>2.4)continue;
                float diJetM = (*thisJet+*thatJet).M();
                if (diJetM>140||diJetM<110) continue;
                HindexList[iList][0] = ij;
                HindexList[iList][1] = jj;
                HindexList[iList][2] = (*thisJet+*thatJet).Pt();
                iList++;
                if (iList>=nCandidates) {
                    cout << "Higgs index out of setting" << endl;
                    return;
                }
                //if (abs(thisJet->Eta()-thatJet->Eta())>1.0) continue;
                //if ((*thisJet+*thatJet).Pt()>HptMax) {
                //}
            }
        } // end of outer loop jet
        sort(HindexList.begin(),HindexList.end(),sortListbyPt);
        if (HindexList[0][0]<0) {h_hNcandi->Fill(0); continue;}
        nPass[1]++;
        
        //cout << "Event: " << jEntry << endl;
        for (int i=0;i<HindexList.size();i++) {
            if (HindexList[i][0]<0) {
                h_hNcandi->Fill(i);
                HindexList.erase(HindexList.begin()+i,HindexList.end());
                break;
            }
            //cout << "jet" << i << " [ " << HindexList[i][0] << " , " << HindexList[i][1] << " , " << HindexList[i][2] << " ]" << endl;
        }
        
        //TLorentzVector* HbJet0 = (TLorentzVector*)genjetP4->At(HindexList[0][0]);
        //TLorentzVector* HbJet1 = (TLorentzVector*)genjetP4->At(HindexList[0][1]);
        //h_hM->Fill((*HbJet0+*HbJet1).M());
        
        // reco A0
        for(int ij=0, iList=0; ij < nGenJet; ij++) {
            TLorentzVector* thisJet = (TLorentzVector*)genjetP4->At(ij);
            if(upeff && thisJet->Pt()<30)continue;
            if(upeff && fabs(thisJet->Eta())>2.4)continue;
            for (int jj=0;jj<ij;jj++) {
                TLorentzVector* thatJet = (TLorentzVector*)genjetP4->At(jj);
                if(upeff && thatJet->Pt()<30) continue;
                if(upeff && fabs(thatJet->Eta())>2.4) continue;
                float diJetM = (*thisJet+*thatJet).M();
                if (diJetM>(A0mass.Atof()+50)||diJetM<(A0mass.Atof()-50)) continue;
                A0indexList[iList][0] = ij;
                A0indexList[iList][1] = jj;
                A0indexList[iList][2] = (*thisJet+*thisJet).Pt();
                iList++;
                if (iList>=nCandidates) {
                    cout << "A0 index out of setting" << endl;
                    return;
                }
                //if (abs(thisJet->Eta()-thatJet->Eta())>1.0) continue;
                //if ((*thisJet+*thatJet).Pt()>A0ptMax) {
                //}
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
            //cout << "jet" << i << " [ " << A0indexList[i][0] << " , " << A0indexList[i][1] << " , " << A0indexList[i][2] << " ]" << endl;
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
                zpM = (*bJet0+*bJet1+*bJet2+*bJet3).M();
                if (zpM>(Zpmass.Atof()+100) || zpM<(Zpmass.Atof()-100)) continue;
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
        if (ZpindexList.size()!=1) continue;
        nPass[4]++;
        
        // apply kinematic selection
        TLorentzVector* HbJet0 = (TLorentzVector*)genjetP4->At(ZpindexList[0][0]);
        TLorentzVector* HbJet1 = (TLorentzVector*)genjetP4->At(ZpindexList[0][1]);
        TLorentzVector* A0bJet0 = (TLorentzVector*)genjetP4->At(ZpindexList[0][2]);
        TLorentzVector* A0bJet1 = (TLorentzVector*)genjetP4->At(ZpindexList[0][3]);
        TLorentzVector* A0recoJet = new TLorentzVector(), *HrecoJet = new TLorentzVector();
        *A0recoJet = *A0bJet0 + *A0bJet1;
        *HrecoJet = *HbJet0 +*HbJet1;
        if ((*HbJet0+*HbJet1).Pt()<250) continue;
        nPass[5]++;
        if ((*A0bJet0+*A0bJet1).Pt()<250) continue;
        nPass[6]++;
        if (A0bJet0->DeltaR(*A0bJet1)>2.2) continue;
        nPass[7]++;
        if (HbJet0->DeltaR(*HbJet1)>1.0) continue;
        nPass[8]++;
        if (abs(HrecoJet->Eta()-A0recoJet->Eta())>0.8) continue;
        nPass[9]++;

        h_hM->Fill((*HbJet0+*HbJet1).M());
        h_hPt->Fill((*HbJet0+*HbJet1).Pt());
        h_hPtAs->Fill(ptAssymetry(HbJet0,HbJet1));
        h_hPtSD->Fill(softDropAs(HbJet0,HbJet1));
        h_hDeltaR->Fill(HbJet0->DeltaR(*HbJet1));
        h_hDeltaPhi->Fill(caldePhi(HbJet0->Phi(),HbJet1->Phi()));
        h_hDeltaEta->Fill(abs(HbJet0->Eta()-HbJet1->Eta()));
        
        h_a0M->Fill((*A0bJet0+*A0bJet1).M());
        h_a0Pt->Fill((*A0bJet0+*A0bJet1).Pt());
        h_a0PtAs->Fill(ptAssymetry(A0bJet0,A0bJet1));
        h_a0PtSD->Fill(softDropAs(A0bJet0,A0bJet1));
        h_a0DeltaR->Fill(A0bJet0->DeltaR(*A0bJet1));
        h_a0DeltaPhi->Fill(caldePhi(A0bJet0->Phi(),A0bJet1->Phi()));
        h_a0DeltaEta->Fill(abs(A0bJet0->Eta()-A0bJet1->Eta()));
        
        h_zpM->Fill((*HbJet0+*HbJet1+*A0bJet0+*A0bJet1).M());
        h_zpPtAs->Fill(ptAssymetry(HrecoJet,A0recoJet));
        h_zpPtSD->Fill(softDropAs(HrecoJet,A0recoJet));
        h_zpDeltaR->Fill(abs(A0recoJet->DeltaR(*HrecoJet)));
        h_zpDeltaPhi->Fill(caldePhi(A0recoJet->Phi(),HrecoJet->Phi()));
        h_zpDeltaEta->Fill(abs(HrecoJet->Eta()-A0recoJet->Eta()));
    } // end of loop over entries
    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    for(int i=0;i<20;i++) if(nPass[i]>0) {
        std::cout << "nPass[" << i << "] = " << nPass[i];
        if (i>0) cout << "\tratio: " << (float)nPass[i]/nPass[i-1] << std::endl;
        else cout << endl; 
    }
    efferr(nPass[4],nTotal);
    efferr(nPass[9],nPass[4]);
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
        h_hNcandi->Draw("hist");
        c1->Print(pdfName.data());
        h_a0Ncandi->Draw("hist");
        c1->Print(pdfName.data());
        h_zpNcandi->Draw("hist");
        c1->Print(pdfNameF.data());
    }
    
    string fileName;
    if (isBG) fileName = Form("QCDbg2HDMbb_%d.root",w);
    else fileName = Form("bb2HDM_MZp%s_MA0%s.root",Zpmass.Data(),A0mass.Data());
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
    h_hNcandi->Write();
    h_a0Ncandi->Write();
    h_zpNcandi->Write();
    outputFile->Close();


}
