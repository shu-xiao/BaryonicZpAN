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

// 10K events TTree --> ~ 600MB space
#define saveTree    false
#define savePDFfile true

using namespace std;
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

struct diJetInfo {
    
    float M;
    float Pt;
    float Assym_pt;
    float Assym_sd;
    diJetInfo(TLorentzVector *v1, TLorentzVector *v2) {
        M = (*v1+*v2).M();
        Pt = (*v1+*v2).Pt();
        Assym_pt = ptAssymetry(v1,v2);
        Assym_sd = softDropAs(v1,v2);
    }
};
struct histList {
    TH1F* h_M;
    TH1F* h_Pt;
    TH1F* h_nComb;
    TH1F* h_ptAs;
    TH1F* h_sdAs;
    TH1F* h_b1Pt;
    TH1F* h_b2Pt;
};
void push_vdiJetInfo(struct histList &sh,vector<struct diJetInfo> &vsdiJet,TLorentzVector *v1,TLorentzVector *v2){
    sh.h_M->Fill((*v1+*v2).M());
    sh.h_Pt->Fill((*v1+*v2).Pt());
    sh.h_ptAs->Fill(ptAssymetry(v1,v2));
    sh.h_sdAs->Fill(softDropAs(v1,v2));
    sh.h_b1Pt->Fill(v1->Pt());
    sh.h_b2Pt->Fill(v2->Pt());
    vsdiJet.push_back(diJetInfo(v1,v2));
}
bool sortByPt(struct diJetInfo s1, struct diJetInfo s2) {return s1.Pt>s2.Pt;}
void sideband(int w=0, std::string inputFile="../../QCDtestBGrootfile/NCUGlobalTuples_76.root") {
//void sideband(int w=0, std::string inputFile="../../QCDtestBGrootfile/NCUGlobalTuples_243.root") {
    

    const float CISVV2CUT_L = 0.5426;
    const float CISVV2CUT_M = 0.8484;
    const float CISVV2CUT_T = 0.9535;

    setNCUStyle(true);
    //gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0001111101.);
   
    vector<struct histList> Ssb_H(3);
    vector<struct histList> Ssb_A0(3);
    vector<vector<struct diJetInfo>> H_diJet(3), A0_diJet(3);  
    //string suffix[3] = {"L","M","T"};
    char *suffix[] = {(char*)"L",(char*)"M",(char*)"T"};

    TFile *f = TFile::Open(Form("Tsideband_%d.root",w),"RECREATE");
    TTree* t = new TTree("tree","tree");
    vector <float> vsb_Hmass_L,  vsb_Hmass_M,  vsb_Hmass_T;
    vector <float> vsb_A0mass_L, vsb_A0mass_M, vsb_A0mass_T;
    if (saveTree) {
        t->Branch("vsb_Hmass_L",  &vsb_Hmass_L);
        t->Branch("vsb_Hmass_M",  &vsb_Hmass_M);
        t->Branch("vsb_Hmass_T",  &vsb_Hmass_T);
        t->Branch("vsb_A0mass_L", &vsb_A0mass_L);
        t->Branch("vsb_A0mass_M", &vsb_A0mass_M);
        t->Branch("vsb_A0mass_T", &vsb_A0mass_T);
    }
    
    const int   nComBinMax_L = 30;
    const int   nComBinMax_MT = 50;
    const float nComMax_L  = nComBinMax_L*10 - 0.5;
    const float nComMax_MT = nComBinMax_MT - 0.5;
    
    TCanvas* c1 = new TCanvas("c1","c1",500,500);
    
    TH1F* h_allEvent   = new TH1F("h_allEvent"  ,"h_allEvent"  , 10, -0.5,9.5);
    for (int i=0;i<3;i++) {
        int nBin = (i==0)?  nComBinMax_L : nComBinMax_MT;
        float max = (i==0)? nComBinMax_L*10-0.5 : nComBinMax_MT-0.5;
        
        Ssb_H[i].h_M        = new TH1F(Form("h_sbHmass_%s",suffix[i]),Form("h_sbHmass_%s",suffix[i]), 100,0,1000); 
        Ssb_A0[i].h_M       = new TH1F(Form("h_sbA0mass_%s",suffix[i]),Form("h_sbA0mass_%s",suffix[i]), 100,0,1000); 
        Ssb_H[i].h_Pt       = new TH1F(Form("h_sbHPt_%s",suffix[i]),Form("h_sbHPt_%s",suffix[i]),100,0,1000);
        Ssb_A0[i].h_Pt      = new TH1F(Form("h_sbA0Pt_%s",suffix[i]),Form("h_sbA0Pt_%s",suffix[i]),100,0,1000);
        Ssb_H[i].h_nComb    = new TH1F(Form("h_nHCom_%s",suffix[i]),Form("h_nHCom_%s",suffix[i]),nBin,-0.5,max); 
        Ssb_A0[i].h_nComb   = new TH1F(Form("h_nA0Com_%s",suffix[i]),Form("h_nA0Com_%s",suffix[i]),nBin,-0.5,max); 
        Ssb_H[i].h_ptAs     = new TH1F(Form("h_H_ptAssymetry_%s",suffix[i]),Form("h_H_ptAssymetry_%s",suffix[i]),100,0,2.5); 
        Ssb_A0[i].h_ptAs    = new TH1F(Form("h_A0_ptAssymetry_%s",suffix[i]),Form("h_A0_ptAssymetry_%s",suffix[i]),100,0,2.5); 
        Ssb_H[i].h_sdAs     = new TH1F(Form("h_H_sdAssymetry_%s",suffix[i]),Form("h_H_sdAssymetry_%s",suffix[i]),60,0,0.6); 
        Ssb_A0[i].h_sdAs    = new TH1F(Form("h_A0_sdAssymetry_%s",suffix[i]),Form("h_A0_sdAssymetry_%s",suffix[i]),60,0,0.6); 
        Ssb_H[i].h_b1Pt     = new TH1F(Form("h_Hb1Pt_%s",suffix[i]),Form("h_Hb1Pt_%s",suffix[i]),100,0,1000);
        Ssb_H[i].h_b2Pt     = new TH1F(Form("h_Hb2Pt_%s",suffix[i]),Form("h_Hb2Pt_%s",suffix[i]),100,0,1000);
        Ssb_A0[i].h_b1Pt    = new TH1F(Form("h_A0b1Pt_%s",suffix[i]),Form("h_A0b1Pt_%s",suffix[i]),100,0,1000);
        Ssb_A0[i].h_b2Pt    = new TH1F(Form("h_A0b2Pt_%s",suffix[i]),Form("h_A0b2Pt_%s",suffix[i]),100,0,1000);
    }
    

    Int_t nPass[20]={0};
    
    TreeReader data(inputFile.data());
    
    int ij, jj, i;
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        
        data.GetEntry(jEntry);
        h_allEvent->Fill(1);
   
        int nGenJet =  data.GetInt("THINnJet");
        if (nGenJet<4) continue;
        for (i=0;i<3;i++) {
            H_diJet[i].clear();
            A0_diJet[i].clear();
        }
        float *vCISVV2 = data.GetPtrFloat("THINjetCISVV2");
        TClonesArray* genjetP4 =  (TClonesArray*) data.GetPtrTObject("THINjetP4");
        vector<bool>& vPassID_L = *(vector<bool>*) data.GetPtr("THINjetPassIDLoose");
        for( ij = 0; ij < nGenJet; ij++ ) {
            
            if (!vPassID_L[ij]) continue;
            
            TLorentzVector* thisJet = (TLorentzVector*)genjetP4->At(ij);
            if ( thisJet->Pt()<30)continue;
            if ( fabs(thisJet->Eta())>2.4)continue;
            
            if ( vCISVV2[ij]<CISVV2CUT_L) continue;
            for ( jj=0 ; jj < ij; jj++ ) {
                
                // veto jet ID
                if (!vPassID_L[jj]) continue;
                
                TLorentzVector* thatJet = (TLorentzVector*)genjetP4->At(jj);
                if ( thatJet->Pt()<30)continue;
                if ( fabs(thatJet->Eta())>2.4)continue;
                
                // cut mass in the sideband region
                static float diJetM = 0;
                diJetM = (*thisJet+*thatJet).M();
                static bool match_H, match_A0;
                match_H  = false;
                match_A0 = false;
                //if (diJetM<=140&&diJetM>=110) continue;
                if (diJetM<=110 || diJetM>=140) match_H = true;
                if (diJetM<=250 || diJetM>=350) match_A0 = true;
                
                // cut on CISVV2 and save into TH1F and pair
                if (match_H  && vCISVV2[ij]>CISVV2CUT_L && vCISVV2[jj]>CISVV2CUT_L) push_vdiJetInfo(Ssb_H[0],  H_diJet[0],  thisJet, thatJet); 
                if (match_H  && vCISVV2[ij]>CISVV2CUT_M && vCISVV2[jj]>CISVV2CUT_M) push_vdiJetInfo(Ssb_H[1],  H_diJet[1],  thisJet, thatJet); 
                if (match_H  && vCISVV2[ij]>CISVV2CUT_T && vCISVV2[jj]>CISVV2CUT_T) push_vdiJetInfo(Ssb_H[2],  H_diJet[2],  thisJet, thatJet);
                if (match_A0 && vCISVV2[ij]>CISVV2CUT_L && vCISVV2[jj]>CISVV2CUT_L) push_vdiJetInfo(Ssb_A0[0], A0_diJet[0], thisJet, thatJet); 
                if (match_A0 && vCISVV2[ij]>CISVV2CUT_M && vCISVV2[jj]>CISVV2CUT_M) push_vdiJetInfo(Ssb_A0[1], A0_diJet[1], thisJet, thatJet); 
                if (match_A0 && vCISVV2[ij]>CISVV2CUT_T && vCISVV2[jj]>CISVV2CUT_T) push_vdiJetInfo(Ssb_A0[2], A0_diJet[2], thisJet, thatJet);
            
            }
    
        }  // end of dijet loop
        
        // fill in nCombination TH1F and sort by Pt
        for (i=0;i<3;i++) {
            if (H_diJet.size()>0)  sort(H_diJet[i].begin()  ,H_diJet[i].end(),  sortByPt);
            if (A0_diJet.size()>0) sort(A0_diJet[i].begin() ,A0_diJet[i].end(), sortByPt);
            Ssb_H[i].h_nComb->Fill(H_diJet[i].size());
            Ssb_A0[i].h_nComb->Fill(A0_diJet[i].size());
        }
        for (i=0;i<H_diJet[0].size();i++)  vsb_Hmass_L.push_back(  (H_diJet[0].size())  ? -99.: H_diJet[0][i].M);
        for (i=0;i<H_diJet[1].size();i++)  vsb_Hmass_M.push_back(  (H_diJet[1].size())  ? -99.: H_diJet[1][i].M);
        for (i=0;i<H_diJet[2].size();i++)  vsb_Hmass_T.push_back(  (H_diJet[2].size())  ? -99.: H_diJet[2][i].M);
        for (i=0;i<A0_diJet[0].size();i++) vsb_A0mass_L.push_back( (A0_diJet[0].size()) ? -99.: A0_diJet[0][i].M);
        for (i=0;i<A0_diJet[1].size();i++) vsb_A0mass_M.push_back( (A0_diJet[1].size()) ? -99.: A0_diJet[1][i].M);
        for (i=0;i<A0_diJet[2].size();i++) vsb_A0mass_T.push_back( (A0_diJet[2].size()) ? -99.: A0_diJet[2][i].M);
        if (saveTree) t->Fill();
         
        // save to root file
    } // end of event loop
    f->Write();
    
    // save as pdf file
    if (savePDFfile) {
        string pdfFileName = Form("Tsideband_QCDbg_%d.pdf",w);
        c1->Print((pdfFileName+"[").data());
        h_allEvent->Draw("hist");
        c1->Print(pdfFileName.data());
        for (i=0;i<3;i++) {
            Ssb_H[i].h_M->Draw("hist");
            c1->Print(pdfFileName.data());
            Ssb_H[i].h_Pt->Draw("hist");
            c1->Print(pdfFileName.data());
            Ssb_H[i].h_nComb->Draw("hist");
            c1->Print(pdfFileName.data());
            Ssb_H[i].h_ptAs->Draw("hist");
            c1->Print(pdfFileName.data());
            Ssb_H[i].h_sdAs->Draw("hist");
            c1->Print(pdfFileName.data());
            Ssb_H[i].h_b1Pt->Draw("hist");
            c1->Print(pdfFileName.data());
            Ssb_H[i].h_b2Pt->Draw("hist");
            c1->Print(pdfFileName.data());
             
        }
        c1->Print((pdfFileName+"]").data());
    }
    f->Close();
    delete f;
}
