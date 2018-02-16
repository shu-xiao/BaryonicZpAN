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
#define savePDFfile false

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

struct diJetInfo(TLorentzVector *v1, TLorentzVector *v2) {
    float M = (*v1+*v2).M();
    float Pt = (*v1+*v2).Pt();
    float Assym_pt = ptAssymetry(v1,v2);
    float Assym_sd = softDropAs(v1,v2);
};
void push_vdiJetInfo(struct diJetInfo &st,TLorentzVector *V1, TLorentzVector *V2) {
    st.vMass.push_back((*V1+*V2).M());
    st.vPt.push_back((*V1+*V2).Pt());
}
void push_vdiJetInfo(vector<pair<float,float>> &vp,TLorentzVector *V1, TLorentzVector *V2) {
    vp.push_back(make_pair((*V1+*V2).M(),(*V1+*V2).Pt()));
}
void push_vdiJetInfo(vector<pair<float,float>> &vp, TH1F* h_M, TLorentzVector *V1, TLorentzVector *V2) {
    vp.push_back(make_pair((*V1+*V2).M(),(*V1+*V2).Pt()));
    h_M->Fill((*V1+*V2).M());
}
bool sortbyPt(pair <float,float> p1, pair <float,float> p2) {return p1.second>p2.second;}
void sideband(int w=0, std::string inputFile="../../QCDtestBGrootfile/NCUGlobalTuples_243.root") {
    

    const float CISVV2CUT_L = 0.5426;
    const float CISVV2CUT_M = 0.8484;
    const float CISVV2CUT_T = 0.9535;

    setNCUStyle(true);
    //gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0001111101.);
   
    typedef vector<pair<float,float>> diJet;
    //struct diJetInfo Ssb_H_L ,Ssb_H_M, Ssb_H_T;
    //struct diJetInfo Ssb_A0_L ,Ssb_A0_M, Ssb_A0_T;
    string suffix[3] = ["L","M","T"];
    diJet Ssb_H_L  ,Ssb_H_M,  Ssb_H_T;
    diJet Ssb_A0_L ,Ssb_A0_M, Ssb_A0_T;
    vector <float> vsb_Hmass_L,  vsb_Hmass_M,  vsb_Hmass_T;
    vector <float> vsb_A0mass_L, vsb_A0mass_M, vsb_A0mass_T;
    
    TFile *f = TFile::Open(Form("Tsideband_%d.root",w),"RECREATE");
    TTree* t = new TTree("tree","tree");
    if (saveTree) {
        t->Branch("vsb_Hmass_L",  &vsb_Hmass_L);
        t->Branch("vsb_Hmass_M",  &vsb_Hmass_M);
        t->Branch("vsb_Hmass_T",  &vsb_Hmass_T);
        t->Branch("vsb_A0mass_L", &vsb_A0mass_L);
        t->Branch("vsb_A0mass_M", &vsb_A0mass_M);
        t->Branch("vsb_A0mass_T", &vsb_A0mass_T);
    }
    
    TCanvas* c1 = new TCanvas("c1","c1",500,500);
    
    TH1F* h_allEvent   = new TH1F("h_allEvent"  ,"h_allEvent"  , 10, -0.5,9.5);
    
    TH1F* h_sbHmass_L  = new TH1F("h_sbHmass_L" ,"h_sbHmass_L" , 100,0,1000);
    TH1F* h_sbHmass_M  = new TH1F("h_sbHmass_M" ,"h_sbHmass_M" , 100,0,1000);
    TH1F* h_sbHmass_T  = new TH1F("h_sbHmass_T" ,"h_sbHmass_T" , 100,0,1000);
    TH1F* h_sbA0mass_L = new TH1F("h_sbA0mass_L","h_sbA0mass_L", 100,0,1000);
    TH1F* h_sbA0mass_M = new TH1F("h_sbA0mass_M","h_sbA0mass_M", 100,0,1000);
    TH1F* h_sbA0mass_T = new TH1F("h_sbA0mass_T","h_sbA0mass_T", 100,0,1000);
    

    const int   nComBinMax_L = 30;
    const int   nComBinMax_MT = 50;
    const float nComMax_L  = nComBinMax_L*10 - 0.5;
    const float nComMax_MT = nComBinMax_MT - 0.5;
    TH1F* h_nHCom_L    = new TH1F("h_nHCom_L"   ,"h_nHCom_L"   , nComBinMax_L ,-0.5,nComMax_L);
    TH1F* h_nHCom_M    = new TH1F("h_nHCom_M"   ,"h_nHCom_M"   , nComBinMax_MT,-0.5,nComMax_MT);
    TH1F* h_nHCom_T    = new TH1F("h_nHCom_T"   ,"h_nHCom_T"   , nComBinMax_MT,-0.5,nComMax_MT);
    TH1F* h_nA0Com_L    = new TH1F("h_nA0Com_L" ,"h_nA0Com_L"  , nComBinMax_L ,-0.5,nComMax_L);
    TH1F* h_nA0Com_M    = new TH1F("h_nA0Com_M" ,"h_nA0Com_M"  , nComBinMax_MT,-0.5,nComMax_MT);
    TH1F* h_nA0Com_T    = new TH1F("h_nA0Com_T" ,"h_nA0Com_T"  , nComBinMax_MT,-0.5,nComMax_MT);

    Int_t nPass[20]={0};
    
    TreeReader data(inputFile.data());
    
    int ij, jj, i;
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        
        data.GetEntry(jEntry);
        h_allEvent->Fill(1);
   
        int nGenJet =  data.GetInt("THINnJet");
        if (nGenJet<4) continue;
        
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
                if (match_H  && vCISVV2[ij]>CISVV2CUT_L && vCISVV2[jj]>CISVV2CUT_L) push_vdiJetInfo(Ssb_H_L,  h_sbHmass_L,  thisJet, thatJet); 
                if (match_H  && vCISVV2[ij]>CISVV2CUT_M && vCISVV2[jj]>CISVV2CUT_M) push_vdiJetInfo(Ssb_H_M,  h_sbHmass_M,  thisJet, thatJet); 
                if (match_H  && vCISVV2[ij]>CISVV2CUT_T && vCISVV2[jj]>CISVV2CUT_T) push_vdiJetInfo(Ssb_H_T,  h_sbHmass_T,  thisJet, thatJet);
                if (match_A0 && vCISVV2[ij]>CISVV2CUT_L && vCISVV2[jj]>CISVV2CUT_L) push_vdiJetInfo(Ssb_A0_L, h_sbA0mass_L, thisJet, thatJet); 
                if (match_A0 && vCISVV2[ij]>CISVV2CUT_M && vCISVV2[jj]>CISVV2CUT_M) push_vdiJetInfo(Ssb_A0_M, h_sbA0mass_M, thisJet, thatJet); 
                if (match_A0 && vCISVV2[ij]>CISVV2CUT_T && vCISVV2[jj]>CISVV2CUT_T) push_vdiJetInfo(Ssb_A0_T, h_sbA0mass_T ,thisJet, thatJet);
            
            }
    
        }  // end of dijet loop
        
        // sort dijet mass vector by pt
        sort(Ssb_H_L.begin(),  Ssb_H_L.end(),  sortbyPt);
        sort(Ssb_H_M.begin(),  Ssb_H_M.end(),  sortbyPt);
        sort(Ssb_H_T.begin(),  Ssb_H_T.end(),  sortbyPt);
        sort(Ssb_A0_L.begin(), Ssb_A0_L.end(), sortbyPt);
        sort(Ssb_A0_M.begin(), Ssb_A0_M.end(), sortbyPt);
        sort(Ssb_A0_T.begin(), Ssb_A0_T.end(), sortbyPt);
        
        // fill in nCombination TH1F
        h_nHCom_L ->Fill(Ssb_H_L.size());
        h_nHCom_M ->Fill(Ssb_H_M.size());
        h_nHCom_T ->Fill(Ssb_H_T.size());
        h_nA0Com_L->Fill(Ssb_A0_L.size());
        h_nA0Com_M->Fill(Ssb_A0_M.size());
        h_nA0Com_T->Fill(Ssb_A0_T.size());

        for (i=0;i<Ssb_H_L.size();i++) vsb_Hmass_L.push_back(Ssb_H_L[i].first);
        for (i=0;i<Ssb_H_M.size();i++) vsb_Hmass_M.push_back(Ssb_H_M[i].first);
        for (i=0;i<Ssb_H_T.size();i++) vsb_Hmass_T.push_back(Ssb_H_T[i].first);
        for (i=0;i<Ssb_A0_L.size();i++) vsb_A0mass_L.push_back(Ssb_A0_L[i].first);
        for (i=0;i<Ssb_A0_M.size();i++) vsb_A0mass_M.push_back(Ssb_A0_M[i].first);
        for (i=0;i<Ssb_A0_T.size();i++) vsb_A0mass_T.push_back(Ssb_A0_T[i].first);
        if (vsb_Hmass_L.size()==0)  vsb_Hmass_L.push_back(-99.);
        if (vsb_Hmass_M.size()==0)  vsb_Hmass_M.push_back(-99.);
        if (vsb_Hmass_T.size()==0)  vsb_Hmass_T.push_back(-99.);
        if (vsb_A0mass_L.size()==0) vsb_A0mass_L.push_back(-99.);
        if (vsb_A0mass_M.size()==0) vsb_A0mass_M.push_back(-99.);
        if (vsb_A0mass_T.size()==0) vsb_A0mass_T.push_back(-99.);
        if (saveTree) t->Fill();
        
        // save to root file
    } // end of event loop
    f->Write();
    
    // save as pdf file
    if (savePDFfile) {
        string pdfFileName = Form("Tsideband_QCDbg_%d.pdf",w);
        c1->Print((pdfFileName+"[").data());
        
        h_sbHmass_L->Draw("hist");
        c1->Print(pdfFileName.data());
        h_nHCom_L->Draw("hist");
        c1->Print(pdfFileName.data());
        h_sbHmass_M->Draw("hist");
        c1->Print(pdfFileName.data());
        h_nHCom_M->Draw("hist");
        c1->Print(pdfFileName.data());
        h_sbHmass_T->Draw("hist");
        c1->Print(pdfFileName.data());
        h_nHCom_T->Draw("hist");
        c1->Print(pdfFileName.data());
        h_sbA0mass_L->Draw("hist");
        c1->Print(pdfFileName.data());
        h_nA0Com_L->Draw("hist");
        c1->Print(pdfFileName.data());
        h_sbA0mass_M->Draw("hist");
        c1->Print(pdfFileName.data());
        h_nA0Com_M->Draw("hist");
        c1->Print(pdfFileName.data());
        h_sbA0mass_T->Draw("hist");
        c1->Print(pdfFileName.data());
        h_nA0Com_T->Draw("hist");
        c1->Print(pdfFileName.data());

        c1->Print((pdfFileName+"]").data());
    }
    f->Close();
    delete f;
}
