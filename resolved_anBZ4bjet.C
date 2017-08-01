
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
using namespace std;
void resolved_anBZ4bjet(int w, std::string inputFile){
    setNCUStyle(true);
    //get TTree from file ...
    TreeReader data(inputFile.data());

    TString zpmass;
    zpmass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*_bb_MZp}; test=${test1%%.root}; echo \"${test}\"",inputFile.data()));

    TCanvas* c1 = new TCanvas("c1","",600,600);

    TH1F* h_passCISt = new TH1F("h_passCISt","h_nJetPassCISSV2_tight",11,-0.5,10.5);
    TH1F* h_passCISm = new TH1F("h_passCISm","h_nJetPassCISSV2_medium",11,-0.5,10.5);
    TH1F* h_passCISl = new TH1F("h_passCISl","h_nJetPassCISSV2_loose",11,-0.5,10.5);
    TH1F* h_hpassCISt0 = new TH1F("h_htobb0_passCISt","h_htobb0_CISSV2_tight",3,-0.5,2.5);
    TH1F* h_hpassCISm0 = new TH1F("h_htobb0_passCISm","h_htobb0_passCISSV2_medium",3,-0.5,2.5);
    TH1F* h_hpassCISl0 = new TH1F("h_htobb0_passCISl","h_htobb0_passCISSV2_loose",3,-0.5,2.5);
    TH1F* h_hpassCISt1 = new TH1F("h_htobb1_passCISt","h_htobb1_passCISSV2_tight",3,-0.5,2.5);
    TH1F* h_hpassCISm1 = new TH1F("h_htobb1_passCISm","h_htobb1_passCISSV2_medium",3,-0.5,2.5);
    TH1F* h_hpassCISl1 = new TH1F("h_htobb1_passCISl","h_htobb1_passCISSV2_loose",3,-0.5,2.5);
    
    TH1F* h_hpt = new TH1F("h_higgsPt", "h_higgsPt_reco", 20,0,1400);


    TH1F* h_oriCIS = new TH1F("h_oriCIS", "h_originCISVV2", 40,0,1);
    TH1F* h_hCIS0 = new TH1F("h_higgsCIS0", "h_higgstobCISVV2_0", 40,0,1);
    TH1F* h_hCIS1 = new TH1F("h_higgsCIS1", "h_higgstobCISVV2_1", 40,0,1);
    
    TH1F* h_hM = new TH1F("h_higgsM", "h_higgsM_reco", 25,0,600);
    TH1F* h_zpM = new TH1F("h_ZpM", "h_ZpM_reco", 25,0,2000);
    TH1F* h_4bM = new TH1F("h_4bM","h_4bM_reco", 25,500,2500);
    
    TH1F* h_hDeltaR = new TH1F("h_HiggstobbdeltaR", "h_HiggstobbDeltaR_reco", 30,0,6);
    TH1F* h_hDeltaPhi = new TH1F("h_HiggstobbdeltaPhi", "h_HiggstobbDeltaPhi_reco", 32,0,3.2);
    TH1F* h_hDeltaEta = new TH1F("h_HiggstobbdeltaEta", "h_HiggstobbDeltaEta_reco", 30,0,6);
    
    Int_t nPass[20]={0};
    

    //for(Long64_t jEntry=0; jEntry<1 ;jEntry++){
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        
        // broken events
        if (jEntry==1370) continue;
        if (jEntry==2974) continue;
        if (jEntry==7594) continue;
        if (jEntry==9888) continue;
        if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        data.GetEntry(jEntry);
        
        //0. has a good vertex
        int nVtx = data.GetInt("nVtx");
        if(nVtx<1)continue;
        nPass[0]++;

        //1. trigger 
        std::string* trigName = data.GetPtrString("hlt_trigName");
        vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));

        bool passTrigger=false;
        for(unsigned int it=0; it< trigResult.size(); it++) {
        
            std::string thisTrig= trigName[it];
            //cout << thisTrig << endl;
            bool results = trigResult[it];
            bool tri[10];
            tri[0] = thisTrig.find("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5")!=std::string::npos;
            tri[1] = thisTrig.find("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50")!=std::string::npos;
            tri[2] = thisTrig.find("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50")!=std::string::npos;
            tri[3] = thisTrig.find("HLT_PFHT800")!=std::string::npos;
            tri[4] = thisTrig.find("HLT_PFHT900")!=std::string::npos;
            tri[5] = thisTrig.find("AK8PFJet360")!=std::string::npos;
            tri[6] = thisTrig.find("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20")!=std::string::npos;
            tri[7] = thisTrig.find("AK8DIPFJET280")!=std::string::npos;
            tri[8] = thisTrig.find("HLT_PFJet260_v9")!=std::string::npos;
            if ( (tri[0] || tri[1] || tri[2] || tri[3] || tri[4] || tri[5] || tri[6] || tri[7] || tri[8] ) && results==1 )    
                {
                    passTrigger=true;
                    break;
                }
        } // end of trigger
        if(!passTrigger) continue;
        nPass[1]++;
        
        // study thin jet CISVV2
        float* thinJetCISV =  data.GetPtrFloat("THINjetCISVV2");
        const int nTHINJets     = data.GetInt("THINnJet");
        TClonesArray* thinjetP4 = (TClonesArray*) data.GetPtrTObject("THINjetP4");
        vector<bool>& passThinJetLooseID = *((vector<bool>*) data.GetPtr("THINjetPassIDLoose"));
        vector<bool>& passThinJetPUID = *((vector<bool>*) data.GetPtr("THINisPUJetIDMedium"));    
    
        int Hindex[2]={-1,-1};
        int Zpindex[2]={-1,-1};
        float maxHpt=-999;
        float maxZppt=-999;
        const float ptCut = 30;
        const float thinCISVV2cutL = 0.5426;
        const float thinCISVV2cutM = 0.8484; // mid cut
        const float thinCISVV2cutT = 0.9535; // tight cut
        // reference: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco 
        int ind = 0;
        int npJetT = 0, npJetM = 0, npJetL = 0;
        float HptMax = -99; 
        for(int ij=0; ij < nTHINJets; ij++) {
            TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(ij);
            if(thisJet->Pt()<ptCut)continue;
            if(fabs(thisJet->Eta())>2.4)continue;
            if(!passThinJetLooseID[ij])continue;
            if(!passThinJetPUID[ij])continue;
            if(thinJetCISV[ij]>thinCISVV2cutT) npJetT++;
            if(thinJetCISV[ij]>thinCISVV2cutM) npJetM++;
            if(thinJetCISV[ij]>thinCISVV2cutL) npJetL++;
            // reco higgs
            for (int jj=0;jj<ij;jj++) {
                TLorentzVector* thatJet = (TLorentzVector*)thinjetP4->At(jj);
                if(thatJet->Pt()<ptCut)continue;
                if(fabs(thatJet->Eta())>2.4)continue;
                if(!passThinJetLooseID[jj])continue;
                if(!passThinJetPUID[jj])continue;
                float diJetM = (*thisJet+*thatJet).M();
                if (diJetM>140||diJetM<110) continue;
                if (abs(thisJet->Eta()-thatJet->Eta())>1.0) continue;
                if (thinJetCISV[jj]>thinCISVV2cutT && thinJetCISV[ij]>thinCISVV2cutT) continue;
                if (thinJetCISV[jj]<thinCISVV2cutL && thinJetCISV[ij]<thinCISVV2cutL) continue;
                if ((*thisJet+*thatJet).Pt()>HptMax) {
                    Hindex[0] = ij;
                    Hindex[1] = jj;
                }
            }
        } // end of outer loop jet
        if (Hindex[0]<0||Hindex[1]<0) continue;
        nPass[2]++;
        TLorentzVector* HbJet0 = (TLorentzVector*)thinjetP4->At(Hindex[0]);
        TLorentzVector* HbJet1 = (TLorentzVector*)thinjetP4->At(Hindex[1]);
        int htobbPl[2]={0}, htobbPm[2]={0}, htobbPt[2]={0};
        for (int i=0;i<2;i++) {
            if (thinJetCISV[Hindex[i]]>thinCISVV2cutL) htobbPl[i]++;
            if (thinJetCISV[Hindex[i]]>thinCISVV2cutM) htobbPm[i]++;
            if (thinJetCISV[Hindex[i]]>thinCISVV2cutT) htobbPt[i]++;
        }
        h_hpassCISt0->Fill(htobbPt[0]);
        h_hpassCISm0->Fill(htobbPm[0]);
        h_hpassCISl0->Fill(htobbPl[0]);
        h_hpassCISt1->Fill(htobbPt[1]);
        h_hpassCISm1->Fill(htobbPm[1]);
        h_hpassCISl1->Fill(htobbPl[1]);
        h_hM->Fill((*HbJet0+*HbJet1).M());
        h_hDeltaR->Fill(HbJet0->DeltaR(*HbJet1));
        h_hDeltaEta->Fill(abs(HbJet0->Eta()-HbJet1->Eta()));
        h_hDeltaPhi->Fill(caldePhi(HbJet0->Phi(),HbJet1->Phi()));
        h_hpt->Fill((*HbJet0+*HbJet1).Pt()); 
        h_passCISt->Fill(npJetT);
        h_passCISm->Fill(npJetM);
        h_passCISl->Fill(npJetL);
    } // end of loop over entries
    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    for(int i=0;i<20;i++) if(nPass[i]>0) std::cout << "nPass[" << i << "] = " << nPass[i] << std::endl;
    efferr(nPass[2],nTotal);
    h_hDeltaR->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf(");
    h_hpt->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_passCISt->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_passCISm->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_passCISl->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf)");
    
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



}
