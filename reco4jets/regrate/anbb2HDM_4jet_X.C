#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1D.h>
#include <TFile.h>
#include "../untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <fstream>
#include <TCanvas.h>
#include "string"
#include "../setNCUStyle.C"
#include "../../gen2HDMsample/genMatch.C"
#include "TMVA_regression_nu_Vali_st.h"
//#include "TMVA_regression_nu_Vali_f.h"

#define saveEffi          true
#define basePtEtaCut    true
#define doGenMatch      false
#define doBTagCut       false
#define CISVV2Cut       0.5426

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
bool sortListbyxSquare(HA0JetInfo s1, HA0JetInfo s2) {return s1.xSqure()<s2.xSqure();}
bool sortListbyWeightxSquare(HA0JetInfo s1, HA0JetInfo s2) {return s1.weightJets()->xSqure()<s2.weightJets()->xSqure();}
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
void savenPass(int nPass[],string fileName) {
    fstream fp;
    fp.open(fileName.data(), ios::out);
    for (int i=0;i<20;i++) {
        if (nPass[i]==0) break;
        fp << (float)nPass[i] << endl;
    }
    fp.close();

}
using namespace std;
void anbb2HDM_4jet_X(int w=0, std::string inputFile="2HDM_MZp1000_MA0300_re.root"){
    
    setNCUStyle(true);
    //gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0001111101.);
    //gStyle->SetOptStat(0);

    //get TTree from file ...
    TreeReader data(inputFile.data());
    vector<vector<int>> genPar = genMatch_base(inputFile.data());
    
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
    //float bin_HT[10] = {50,100,200,300,500,700,1000,1500,2000,3000};
    //TH1F* h_HT = new TH1F("h_HT","h_HT",9,bin_HT);
    TH1F* h_HT = new TH1F("h_HT","h_HT",60,0,3000);
    
    TH1F* h_hPtAs = new TH1F("h_hPtAs","h_higgsPtAssymetry",40,0,2);
    TH1F* h_hPtAs_ori = new TH1F("h_hPtAs_ori","h_higgsPtAssymetry",40,0,2);
    TH1F* h_a0PtAs = new TH1F("h_a0PtAs","h_A0PtAssymetry",40,0,2);
    TH1F* h_zpPtAs = new TH1F("h_zpPtAs","h_ZpPtAssymetry",60,0,3.0);
    TH1F* h_zpPtAs_ori = new TH1F("h_zpPtAs_ori","h_ZpPtAssymetry",60,0,3.0);

    TH1F* h_hPtSD = new TH1F("h_hPtSD","h_higgsPtSD",40,0,0.7);
    TH1F* h_hPtSD_ori = new TH1F("h_hPtSD_ori","h_higgsPtSD",40,0,0.7);
    TH1F* h_a0PtSD = new TH1F("h_a0PtSD","h_A0PtSD",30,0,0.6);
    TH1F* h_zpPtSD = new TH1F("h_zpPtSD","h_ZpPtSD",30,0.1,0.7);
    TH1F* h_zpPtSD_ori = new TH1F("h_zpPtSD_ori","h_ZpPtSD",30,0.1,0.7);
    
    const int njets = 20;
    const int ncom = 50;
    TH1F* h_NgoodJets = new TH1F("h_NgoodJets", "h_NgoodJets", njets,-0.5,njets-0.5);
    TH1F* h_Ncom = new TH1F("h_Ncom", "h_Ncom", ncom,-0.5,ncom-0.5);

    TH1F* h_hM = new TH1F("h_higgsM", "h_higgsM_weight", 40,0,200);
    TH1F* h_hM_ori = new TH1F("h_higgsM_ori", "h_higgsM_reco", 40,0,200);
    TH1F* h_hM_all = new TH1F("h_higgsM_all", "h_higgsM_all", 100,50,150);
    //TH1F* h_a0M = new TH1F("h_a0M", "h_A0M_weight", 24,A0mass.Atof()-60,A0mass.Atof()+60);
    TH1F* h_a0M = new TH1F("h_a0M", "h_A0M_weight", 100,0,500);
    TH1F* h_a0M_ori = new TH1F("h_A0M_ori", "h_A0M_ori", 100,0,500);
    TH1F* h_zpM = new TH1F("h_ZpM", "h_ZpM_weight", 24,Zpmass.Atof()-120,Zpmass.Atof()+120);
    TH1F* h_zpM_ori = new TH1F("h_ZpM_ori", "h_ZpM_reco", 24,Zpmass.Atof()-120,Zpmass.Atof()+120);
    
    TH1F* h_hPt = new TH1F("h_higgsPt", "h_higgsPt_reco", 60,0,1200);
    TH1F* h_hPt_ori = new TH1F("h_higgsPt_ori", "h_higgsPt_reco", 60,0,1200);
    TH1F* h_a0Pt = new TH1F("h_a0Pt", "h_A0Pt_reco", 60,0,1200);
    
    TH1F* h_hDeltaR = new TH1F("h_HiggstobbDeltaR", "h_HiggstobbDeltaR_reco", 40,0,4);
    TH1F* h_hDeltaR_ori = new TH1F("h_HiggstobbDeltaR_ori", "h_HiggstobbDeltaR_reco", 40,0,4);
    TH1F* h_a0DeltaR = new TH1F("h_A0tobbDeltaR", "h_A0tobbDeltaR_reco", 40,0,4);
    TH1F* h_zpDeltaR = new TH1F("h_ZptoHA0DeltaR", "h_ZptoHA0DeltaR_reco", 35,1,4.5);
    TH1F* h_zpDeltaR_ori = new TH1F("h_ZptoHA0DeltaR_ori", "h_ZptoHA0DeltaR_reco", 35,1,4.5);
    TH1F* h_hDeltaEta = new TH1F("h_HiggstobbDeltaEta", "h_HiggstobbDeltaEta_reco", 40,0,4);
    TH1F* h_hDeltaEta_ori = new TH1F("h_HiggstobbDeltaEta_ori", "h_HiggstobbDeltaEta_reco", 40,0,4);
    TH1F* h_a0DeltaEta = new TH1F("h_A0tobbDeltaEta", "h_A0tobbDeltaEta_reco", 40,0,4);
    TH1F* h_zpDeltaEta_ori = new TH1F("h_ZptoHA0DeltaEta_ori", "h_ZptoHA0DeltaEta_reco", 40,0,4);
    TH1F* h_zpDeltaEta = new TH1F("h_ZptoHA0DeltaEta", "h_ZptoHA0DeltaEta_reco", 40,0,4);
    TH1F* h_hDeltaPhi = new TH1F("h_HiggstobbDeltaPhi", "h_HiggstobbDeltaPhi_reco", 32,0,3.2);
    TH1F* h_hDeltaPhi_ori = new TH1F("h_HiggstobbDeltaPhi_ori", "h_HiggstobbDeltaPhi_reco", 32,0,3.2);
    TH1F* h_a0DeltaPhi = new TH1F("h_A0tobbDeltaPhi", "h_A0obbDeltaPhi_reco", 32,0,3.2);
    TH1F* h_zpDeltaPhi = new TH1F("h_ZptoHA0DeltaPhi", "h_ZptoHA0DeltaPhi_reco", 32,0,3.2);
    TH1F* h_zpDeltaPhi_ori = new TH1F("h_ZptoHA0DeltaPhi_ori", "h_ZptoHA0DeltaPhi_reco", 32,0,3.2);
    
    TH1F *h_CISVV2_H[2], *h_CISVV2_A0[2];
    TH1F *h_TMVAweight_sw[2];
    TH1F *h_TMVAweight_ws[2];
    h_hM->SetLineColor(kBlue);
    h_zpM->SetLineColor(kBlue);
    h_hPt->SetLineColor(kBlue);
    h_hDeltaR->SetLineColor(kBlue);
    h_hDeltaPhi->SetLineColor(kBlue);
    h_hDeltaEta->SetLineColor(kBlue);
    h_hPtAs->SetLineColor(kBlue);
    h_hPtSD->SetLineColor(kBlue);
    h_zpDeltaR->SetLineColor(kBlue);
    h_zpDeltaPhi->SetLineColor(kBlue);
    h_zpDeltaR->SetLineColor(kBlue);
    h_zpDeltaEta->SetLineColor(kBlue);
    h_zpPtSD->SetLineColor(kBlue);
    h_zpPtAs->SetLineColor(kBlue);
    for (int i=0;i<2;i++) {
        h_CISVV2_H[i]= new TH1F(Form("h_CISVV2_Htobb_%d",i),Form("h_CISVV2_Htobb_%d",i),22,0,1.1);
        h_CISVV2_A0[i]= new TH1F(Form("h_CISVV2_A0tobb_%d",i),Form("h_CISVV2_A0tobb_%d",i),22,0,1.1);
        h_TMVAweight_ws[i] = new TH1F(Form("TMVAweightWS_%d",i),Form("TMVAweightWS_%d",i),40,0.0,2);
        h_TMVAweight_sw[i] = new TH1F(Form("TMVAweightSW_%d",i),Form("TMVAweightSW_%d",i),40,0.0,2);
        
    }
    Int_t nPass[20]={0};
    int maxHIndexNum = 0, maxZpIndexNum = 0;
    
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        
        if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        data.GetEntry(jEntry);
        //if (jEntry>50) break; 
        nPass[0]++;
        float HT = data.GetFloat("HT");
        h_HT->Fill(HT);
        h_allEvent->Fill(1);
        
        //0. has a good vertex
        int nGenJet =  data.GetInt("THINnJet");
        if(nGenJet<4) continue;
        nPass[1]++;
        
        vector<bool>& vPassID_L = *(vector<bool>*) data.GetPtr("THINjetPassIDLoose");
        TClonesArray* genjetP4 =  (TClonesArray*) data.GetPtrTObject("THINjetP4");
        float *CISVV2 = data.GetPtrFloat("THINjetCISVV2");
        
        vector<vector<float>> HindexList, A0indexList, ZpindexList;
        vector<TLorentzVector*> genHA0Par;
        
        if (!isBG||doGenMatch) {
            TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
            for (int i=0;i<4;i++) genHA0Par.push_back((TLorentzVector*)genParP4->At(genPar[jEntry][i]));
        }
        
        // find good jet
        static vector<int> goodJet;
        goodJet.clear();
        static TLorentzVector* thisJet, *thatJet;
        for (int i=0;i<nGenJet;i++) {
            if (!vPassID_L[i]) continue;
            if (CISVV2[i]<CISVV2Cut) continue;
            thisJet = (TLorentzVector*)genjetP4->At(i);
            if (thisJet->Pt()<30) continue;
            if (thisJet->Eta()>2.4) continue;
            goodJet.push_back(i);
        }
        h_NgoodJets->Fill(goodJet.size());
        if (goodJet.size()<4) continue;
        nPass[2]++;
        TMVAinputInfo info(data);
        /*
        HA0JetInfo    fourJet(info,goodJet[0],goodJet[1],goodJet[2],goodJet[3],1);
        thisJet = (TLorentzVector*)genjetP4->At(goodJet[0]);
        thatJet = (TLorentzVector*)genjetP4->At(goodJet[1]);
        float weight[2];
        weight[0] = TMVA_15plus3_jetGenJet_nu_3_1(info,goodJet[0],thisJet->DeltaR(*thatJet));
        weight[1] = TMVA_15plus3_jetGenJet_nu_3_1(info,goodJet[1],thisJet->DeltaR(*thatJet));
        cout << weight[0] << " " << fourJet.getWeight()[0] <<"\t" << weight[1] << " " << fourJet.getWeight()[1] << endl;
        */
        
        static vector<HA0JetInfo> fourJetList;
        static HA0JetInfo minX2, minWeightX2; 
        fourJetList.clear();
        static int hi, hj, ai, aj;
        cout << "good" << endl;
        for (hi=0;hi<goodJet.size();hi++) {
            for (hj=0;hj<hi;hj++) {
                for (ai=0;ai<goodJet.size();ai++) {
                    if ( ai==hi || ai==hj ) continue;
                    for (aj=0;aj<ai;aj++) {
                        if ( aj==hi || aj==hj ) continue;
                        fourJetList.push_back(HA0JetInfo(info,hi,hj,ai,aj,1));
                    }
                }
            }
        }
        cout << "loop" << endl;
        h_Ncom->Fill(fourJetList.size());
        if (fourJetList.size()==0) continue;
        nPass[3]++;

        sort(fourJetList.begin(),fourJetList.end(),sortListbyxSquare);
        minX2       = fourJetList[0];
        sort(fourJetList.begin(),fourJetList.end(),sortListbyWeightxSquare);
        minWeightX2 = fourJetList[0];
        
        h_hM_ori->Fill(minX2.Mh());
        h_hM->Fill(minWeightX2.weightJets()->Mh());
        
        h_a0M->Fill(minX2.MA0());
        h_a0M_ori->Fill(minWeightX2.weightJets()->MA0());
        
        h_zpM->Fill(minX2.MZp());
        h_zpM_ori->Fill(minWeightX2.weightJets()->MZp());
        h_TMVAweight_sw[0]->Fill(minX2.getWeight(0)); 
        h_TMVAweight_sw[1]->Fill(minX2.getWeight(1)); 
        h_TMVAweight_ws[0]->Fill(minWeightX2.getWeight(0)); 
        h_TMVAweight_ws[1]->Fill(minWeightX2.getWeight(1));
        //h_TMVAweight_ws[0]->Fill(minX2.weightJets()->getWeight(0)); 
        //h_TMVAweight_ws[1]->Fill(minX2.weightJets()->getWeight(1));
        //cout << minX2.getWeight(0) << " " << minX2.getWeight(1)  << "| \t" << minX2.weightJets()->getWeight(0) << " " << minX2.weightJets()->getWeight(1) << endl;
        //cout << minX2.getWeight(0) << " " << minX2.getWeight(1)  << "| \t" << minWeightX2.getWeight(0) << " " << minWeightX2.getWeight(1) << endl;
        
        
        
        
        
        
        
        
        
        
       /* 
        TLorentzVector* HbJet0_base = (TLorentzVector*)genjetP4->At(ZpindexList[0][0]);
        TLorentzVector* HbJet1_base = (TLorentzVector*)genjetP4->At(ZpindexList[0][1]);
        static float TMVAweight[2];
        static float TMVAweight_tt[2];
        TMVAweight[0] = TMVA_15plus3_jetGenJet_nu_3_1(data,(int)ZpindexList[0][0],HbJet0_base->DeltaR(*HbJet1_base));
        TMVAweight[1] = TMVA_15plus3_jetGenJet_nu_3_1(data,(int)ZpindexList[0][1],HbJet0_base->DeltaR(*HbJet1_base));
        TMVAweight_tt[0] = TMVA_15plus3_jetGenJet_nu_3_1(info,(int)ZpindexList[0][0],HbJet0_base->DeltaR(*HbJet1_base));
        TMVAweight_tt[1] = TMVA_15plus3_jetGenJet_nu_3_1(info,(int)ZpindexList[0][1],HbJet0_base->DeltaR(*HbJet1_base));
        cout << TMVAweight[0] << " " << TMVAweight_tt[0] << "\t" << TMVAweight[1] << " " << TMVAweight_tt[1] << endl;
        static TLorentzVector v1 = TLorentzVector(),v2 = TLorentzVector();
        v1.SetPtEtaPhiE(HbJet0_base->Pt()*TMVAweight[0],HbJet0_base->Eta(),HbJet0_base->Phi(),HbJet0_base->E()*TMVAweight[0]); 
        v2.SetPtEtaPhiE(HbJet1_base->Pt()*TMVAweight[1],HbJet1_base->Eta(),HbJet1_base->Phi(),HbJet1_base->E()*TMVAweight[1]); 
        TLorentzVector* HbJet0 = &v1;
        TLorentzVector* HbJet1 = &v2;
        
        
        // Draw figures
        h_hM->Fill((*HbJet0+*HbJet1).M());
        h_hPt->Fill((*HbJet0+*HbJet1).Pt());
        h_hPtAs->Fill(ptAssymetry(HbJet0,HbJet1));
        h_hPtSD->Fill(softDropAs(HbJet0,HbJet1));
        h_hDeltaR->Fill(HbJet0->DeltaR(*HbJet1));
        h_hDeltaPhi->Fill(caldePhi(HbJet0->Phi(),HbJet1->Phi()));
        h_hDeltaEta->Fill(abs(HbJet0->Eta()-HbJet1->Eta()));

        h_hM_ori->Fill((*HbJet0_base+*HbJet1_base).M());
        h_hPt_ori->Fill((*HbJet0_base+*HbJet1_base).Pt());
        h_hPtSD_ori->Fill(softDropAs(HbJet0_base,HbJet1_base));
        h_hPtAs_ori->Fill(ptAssymetry(HbJet0_base,HbJet1_base));
        h_hDeltaR_ori->Fill(HbJet0_base->DeltaR(*HbJet1_base));
        h_hDeltaPhi_ori->Fill(caldePhi(HbJet0_base->Phi(),HbJet1_base->Phi()));
        h_hDeltaEta_ori->Fill(abs(HbJet0_base->Eta()-HbJet1_base->Eta()));
        
        h_TMVAweight[0]->Fill(TMVAweight[0]);
        h_TMVAweight[1]->Fill(TMVAweight[1]);

        TLorentzVector* A0bJet0 = (TLorentzVector*)genjetP4->At(ZpindexList[0][2]);
        TLorentzVector* A0bJet1 = (TLorentzVector*)genjetP4->At(ZpindexList[0][3]);
        h_a0M->Fill((*A0bJet0+*A0bJet1).M());
        h_a0Pt->Fill((*A0bJet0+*A0bJet1).Pt());
        h_a0PtAs->Fill(ptAssymetry(A0bJet0,A0bJet1));
        h_a0PtSD->Fill(softDropAs(A0bJet0,A0bJet1));
        h_a0DeltaR->Fill(A0bJet0->DeltaR(*A0bJet1));
        h_a0DeltaPhi->Fill(caldePhi(A0bJet0->Phi(),A0bJet1->Phi()));
        h_a0DeltaEta->Fill(abs(A0bJet0->Eta()-A0bJet1->Eta()));
        
        TLorentzVector* A0recoJet = new TLorentzVector(), *HrecoJet = new TLorentzVector();
        TLorentzVector*HrecoJet_base = new TLorentzVector();
        *A0recoJet = *A0bJet0 + *A0bJet1;
        *HrecoJet = *HbJet0 +*HbJet1;
        *HrecoJet_base = *HbJet0_base +*HbJet1_base;
        
        h_zpM_ori->Fill((*HbJet0_base+*HbJet1_base+*A0bJet0+*A0bJet1).M());
        h_zpPtAs_ori->Fill(ptAssymetry(HrecoJet_base,A0recoJet));
        h_zpPtSD_ori->Fill(softDropAs(HrecoJet_base,A0recoJet));
        h_zpDeltaR_ori->Fill(abs(A0recoJet->DeltaR(*HrecoJet_base)));
        h_zpDeltaPhi_ori->Fill(caldePhi(A0recoJet->Phi(),HrecoJet_base->Phi()));
        h_zpDeltaEta_ori->Fill(abs(HrecoJet_base->Eta()-A0recoJet->Eta()));
        
        h_zpM->Fill((*HbJet0+*HbJet1+*A0bJet0+*A0bJet1).M());
        h_zpPtAs->Fill(ptAssymetry(HrecoJet,A0recoJet));
        h_zpPtSD->Fill(softDropAs(HrecoJet,A0recoJet));
        h_zpDeltaR->Fill(abs(A0recoJet->DeltaR(*HrecoJet)));
        h_zpDeltaPhi->Fill(caldePhi(A0recoJet->Phi(),HrecoJet->Phi()));
        h_zpDeltaEta->Fill(abs(HrecoJet->Eta()-A0recoJet->Eta()));
    

        delete A0recoJet;
        delete HrecoJet;
        delete HrecoJet_base; 
        */
    } // end of loop over entries
    /*
    h_zpPtSD->SetMaximum(1000);
    h_zpPtAs_ori->SetMaximum(150); 
    */
    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    
    for(int i=0;i<20;i++) if(nPass[i]>0) std::cout << "nPass[" << i << "] = " << nPass[i] << std::endl;
    efferr(nPass[3],nTotal);
    TLegend* leg = new TLegend(0.68,0.8,0.9,0.9);
    leg->AddEntry(h_hM,"TMVA Reweight");
    leg->AddEntry(h_hM_ori,"origin");
    if (!isBG) 
    {
        string pdfName;
        if (doGenMatch) pdfName = Form("anGenMatchJet_bb2HDM_MZp%s_MA0%s.pdf",Zpmass.Data(),A0mass.Data());
        else if (doBTagCut) pdfName = Form("anRecoBTagJet_bb2HDM_MZp%s_MA0%s.pdf",Zpmass.Data(),A0mass.Data());
        else pdfName = Form("anRecoJet_bb2HDM_MZp%s_MA0%s.pdf",Zpmass.Data(),A0mass.Data());
        c1->Print((pdfName+"[").data());
        h_HT->Draw("hist");
        c1->Print(pdfName.data());
        h_NgoodJets->Draw("hist");
        c1->Print(pdfName.data());
        h_Ncom->Draw("hist");
        c1->Print(pdfName.data());
        h_zpM_ori->Draw("hist");
        h_zpM->Draw("histsame");
        leg->Draw();
        c1->Print(pdfName.data());
        h_hM_ori->Draw("hist");
        h_hM->Draw("histsame");
        leg->Draw();
        c1->Print(pdfName.data());
        h_a0M->Draw("hist");
        h_a0M_ori->Draw("hist");
        c1->Print(pdfName.data());
        /*
        h_hPt->Draw("hist");
        h_hPt_ori->Draw("histsame");
        leg->Draw();
        c1->Print(pdfName.data());
        h_a0Pt->Draw("hist");
        c1->Print(pdfName.data());
        h_hPtAs->Draw("hist");
        h_hPtAs_ori->Draw("histsame");
        leg->Draw();
        c1->Print(pdfName.data());
        h_a0PtAs->Draw("hist");
        c1->Print(pdfName.data());
        h_zpPtAs_ori->Draw("hist");
        h_zpPtAs->Draw("histsame");
        leg->Draw();
        c1->Print(pdfName.data());
        h_hPtSD->Draw("hist");
        h_hPtSD_ori->Draw("histsame");
        leg->Draw();
        c1->Print(pdfName.data());
        h_a0PtSD->Draw("hist");
        c1->Print(pdfName.data());
        h_zpPtSD->Draw("hist");
        h_zpPtSD_ori->Draw("histsame");
        leg->Draw();
        c1->Print(pdfName.data());
        //h_hDeltaR_ori->Draw("hist");
        h_hDeltaR->Draw("hist");
        h_hDeltaR_ori->Draw("histsame");
        leg->Draw();
        c1->Print(pdfName.data());
        h_a0DeltaR->Draw("hist");
        c1->Print(pdfName.data());
        h_zpDeltaR->Draw("hist");
        h_zpDeltaR_ori->Draw("histsame");
        leg->Draw();
        c1->Print(pdfName.data());
        h_hDeltaEta->Draw("hist");
        h_hDeltaEta_ori->Draw("histsame");
        leg->Draw();
        c1->Print(pdfName.data());
        h_a0DeltaEta->Draw("hist");
        c1->Print(pdfName.data());
        h_zpDeltaEta->Draw("hist");
        h_zpDeltaEta_ori->Draw("histsame");
        leg->Draw();
        c1->Print(pdfName.data());
        h_hDeltaPhi->Draw("hist");
        h_hDeltaPhi_ori->Draw("histsame");
        leg->Draw();
        c1->Print(pdfName.data());
        h_a0DeltaPhi->Draw("hist");
        c1->Print(pdfName.data());
        h_zpDeltaPhi->Draw("hist");
        h_zpDeltaPhi_ori->Draw("histsame");
        leg->Draw();
        c1->Print(pdfName.data());
        h_CISVV2_H[0]->Draw("hist");
        c1->Print(pdfName.data());
        h_CISVV2_H[1]->Draw("hist");
        c1->Print(pdfName.data());
        h_CISVV2_A0[0]->Draw("hist");
        c1->Print(pdfName.data());
        h_CISVV2_A0[1]->Draw("hist");
        c1->Print(pdfName.data());
        h_hNcandi->Draw("hist text0");
        c1->Print(pdfName.data());
        h_a0Ncandi->Draw("hist text0");
        c1->Print(pdfName.data());
        h_zpNcandi->Draw("hist text0");
        c1->Print(pdfName.data());
        */
        h_TMVAweight_sw[0]->Draw("hist");
        h_TMVAweight_sw[1]->Draw("histsame");
        c1->Print(pdfName.data());
        h_TMVAweight_ws[0]->Draw("hist");
        h_TMVAweight_ws[1]->Draw("histsame");
        c1->Print(pdfName.data());
        c1->Print((pdfName+"]").data());
    }
    /*
    string fileName; 
    if (isBG) fileName = Form("QCDbg2HDMbb_%d.root",w);
    else if (doBTagCut) fileName = Form("sigRootFile/bb2HDM_recoBTag_MZp%s_MA0%s.root",Zpmass.Data(),A0mass.Data());
    else fileName = Form("sigRootFile/bb2HDM_reco_MZp%s_MA0%s.root",Zpmass.Data(),A0mass.Data());
    if (!saveEffi&&false) {
        TFile* outputFile = new TFile(fileName.data(),"recreate");
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
        h_CISVV2_H[0]->Write();
        h_CISVV2_H[1]->Write();
        h_CISVV2_A0[0]->Write();
        h_CISVV2_A0[1]->Write();
        h_hNcandi->Write();
        h_a0Ncandi->Write();
        h_zpNcandi->Write();
        
        outputFile->Close();
    }
    string fileNPassName;
    if (doBTagCut) fileNPassName = Form("../njetsnPass/effi_Zpmass%s_A0mass%s_recoBTag4jets.txt",Zpmass.Data(),A0mass.Data());
    else fileNPassName = Form("../njetsnPass/effi_Zpmass%s_A0mass%s_reco4jets.txt",Zpmass.Data(),A0mass.Data());
    savenPass(nPass,fileNPassName);
    if (saveEffi) {
        string effifilename;
        if (doBTagCut) effifilename = Form("../njetsEffi/effi_Zpmass%s_A0mass%s_recoBTag4jets.txt",Zpmass.Data(),A0mass.Data());
        else effifilename = Form("../njetsEffi/effi_Zpmass%s_A0mass%s_reco4jets.txt",Zpmass.Data(),A0mass.Data());
        //string nPassfilename = Form("../njetsnPass/effi_Zpmass%s_A0mass%s_4jets.txt",Zpmass.Data(),A0mass.Data());
        fstream fp;
        fp.open(effifilename.data(), ios::out);
        fp << (float)nPass[4]/nTotal << endl;
        fp.close();
        
    }*/
}
