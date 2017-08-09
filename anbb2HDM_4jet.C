
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

using namespace std;
void anbb2HDM_4jet(int w, std::string inputFile){
    setNCUStyle(true);
    //get TTree from file ...
    TreeReader data(inputFile.data());
    
    bool upeff = true;
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
    
    TH1F* h_hNcandi = new TH1F("h_hNcandi", "h_higgs_NcandidatePairJets", 10,-0.5,9.5);
    TH1F* h_a0Ncandi = new TH1F("h_a0Ncandi", "h_A0_NcandidatePairJets", 10,-0.5,9.5);
    TH1F* h_zpNcandi = new TH1F("h_zpNcandi", "h_Zp_NcandidatePairJets", 10,-0.5,9.5);

    TH1F* h_hM = new TH1F("h_higgsM", "h_higgsM_genJet", 50,100,150);
    TH1F* h_a0M = new TH1F("h_a0M", "h_A0M_genJet", 24,240,360);
    TH1F* h_zpM = new TH1F("h_ZpM", "h_ZpM_genJet", 24,1380,1620);
    
    TH1F* h_hDeltaR = new TH1F("h_HiggstobbdeltaR", "h_HiggstobbDeltaR_reco", 30,0,6);
    TH1F* h_hDeltaPhi = new TH1F("h_HiggstobbdeltaPhi", "h_HiggstobbDeltaPhi_reco", 32,0,3.2);
    TH1F* h_hDeltaEta = new TH1F("h_HiggstobbdeltaEta", "h_HiggstobbDeltaEta_reco", 30,0,6);
    
    Int_t nPass[20]={0};
    
    //for(Long64_t jEntry=0; jEntry<100 ;jEntry++){
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        
        // broken events
        if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        data.GetEntry(jEntry);
        
        //0. has a good vertex
        int nGenJet = data.GetInt("ak4nGenJet");
        if(nGenJet<4) continue;
        nPass[0]++;
        
        TClonesArray* genjetP4 = (TClonesArray*) data.GetPtrTObject("ak4GenJetP4");
        int* genMomParId= data.GetPtrInt("genMomParId");
        float HptMax = -99, A0ptMax = -99; 
        vector<vector<float>> HindexList, A0indexList, ZpindexList;
        vector<float> row(5,-99);
        HindexList.assign(30,row);
        A0indexList.assign(30,row);
        ZpindexList.assign(30,row);
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
                if (diJetM>145||diJetM<105) continue;
                HindexList[iList][0] = ij;
                HindexList[iList][1] = jj;
                HindexList[iList][2] = (*thisJet+*thatJet).Pt();
                if (iList>=29) {
                    cout << "Higgs index out of setting" << endl;
                    return;
                }
                iList++;
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
                if (diJetM>350||diJetM<250) continue;
                A0indexList[iList][0] = ij;
                A0indexList[iList][1] = jj;
                A0indexList[iList][2] = (*thisJet+*thisJet).Pt();
                if (iList>=29) {
                    cout << "A0 index out of setting" << endl;
                    return;
                }
                iList++;
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
                zpM = (*bJet0+*bJet1+*bJet2+*bJet3).M();
                if (zpM>1600 || zpM<1400) continue;
                ZpindexList[iList][0] = HindexList[i][0];
                ZpindexList[iList][1] = HindexList[i][1];
                ZpindexList[iList][2] = A0indexList[j][0];
                ZpindexList[iList][3] = A0indexList[j][1];
                ZpindexList[iList][5] = (*bJet0+*bJet1+*bJet2+*bJet3).Pt();
                if (iList>=29) {
                    cout << "Zp index out of setting" << endl;
                    return;
                }
                iList++;
            }
        }
        if (ZpindexList[0][0]<0) {h_zpNcandi->Fill(0); continue;}
        nPass[3]++;
        sort(ZpindexList.begin(),ZpindexList.end(),sortListbyPt);
        for (int i=0;i<ZpindexList.size();i++) {
            if (ZpindexList[i][0]<0) {
                h_zpNcandi->Fill(i);
                ZpindexList.erase(ZpindexList.begin()+i,ZpindexList.end());
                break;
            }
            //cout << "jet" << i << " [ " << ZpindexList[i][0] << " , " << ZpindexList[i][1] << " , " << ZpindexList[i][2] << " ]" << endl;
        }
        
        TLorentzVector* HbJet0 = (TLorentzVector*)genjetP4->At(ZpindexList[0][0]);
        TLorentzVector* HbJet1 = (TLorentzVector*)genjetP4->At(ZpindexList[0][1]);
        h_hM->Fill((*HbJet0+*HbJet1).M());
        TLorentzVector* A0bJet0 = (TLorentzVector*)genjetP4->At(ZpindexList[0][2]);
        TLorentzVector* A0bJet1 = (TLorentzVector*)genjetP4->At(ZpindexList[0][3]);
        h_a0M->Fill((*A0bJet0+*A0bJet1).M());
        h_zpM->Fill((*HbJet0+*HbJet1+*A0bJet0+*A0bJet1).M());
    } // end of loop over entries
    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    for(int i=0;i<20;i++) if(nPass[i]>0) std::cout << "nPass[" << i << "] = " << nPass[i] << std::endl;
    efferr(nPass[3],nTotal);
    h_zpM->Draw("hist");
    c1->Print("anGenJet_bb2HDM_MZp1500_MA0300.pdf(");
    h_hM->Draw("hist");
    c1->Print("anGenJet_bb2HDM_MZp1500_MA0300.pdf");
    h_a0M->Draw("hist");
    c1->Print("anGenJet_bb2HDM_MZp1500_MA0300.pdf");
    h_hNcandi->Draw("hist");
    c1->Print("anGenJet_bb2HDM_MZp1500_MA0300.pdf");
    h_a0Ncandi->Draw("hist");
    c1->Print("anGenJet_bb2HDM_MZp1500_MA0300.pdf");
    h_zpNcandi->Draw("hist");
    c1->Print("anGenJet_bb2HDM_MZp1500_MA0300.pdf)");
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
