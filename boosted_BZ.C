
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
using namespace std;
void boosted_BZ(int w, std::string inputFile){
    setNCUStyle(true);
    //get TTree from file ...
    TreeReader data(inputFile.data());

    TString zpmass;
    zpmass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*_bb_MZp}; test=${test1%%.root}; echo \"${test}\"",inputFile.data()));
    string outputFile=Form("effiMZp%s.txt",zpmass.Data());

    TCanvas* c1 = new TCanvas("c1","",889*1.5,768);
    TH1F* h_4bJetM = new TH1F("h_4bMass", "4b mass", 100,100,3300);
    TH1F* h_HT = new TH1F("h_HT", "HT", 20,0,1500);
    TH1F* h_jetM[2], *h_jetPt[2], *h_jetEta[2];
    for (int i=0;i<2;i++) {
        h_jetM[i] = new TH1F(Form("h_jetM_%d",i),Form("h_jetM_%d",i),20,0,1000);
        h_jetPt[i] = new TH1F(Form("h_jetPt_%d",i),Form("h_jetPt_%d",i),20,0,1000);
        h_jetEta[i] = new TH1F(Form("h_jetEta_%d",i),Form("h_jetEta_%d",i),20,-3,3);
    }

    TH1F* h_ptb = new TH1F("h_ptb0", "h_ptb0", 20,0,1000);
    TH1F* h_pta = new TH1F("h_pta0", "h_pta0", 20,0,1000);
    TH1F* h_ptb1 = new TH1F("h_ptb1", "h_ptb1", 20,0,1000);
    TH1F* h_pta1 = new TH1F("h_pta1", "h_pta1", 20,0,1000);

    TH1F* h_taub = new TH1F("h_taub0", "h_taub0", 20,0,1);
    TH1F* h_taua = new TH1F("h_taua0", "h_taua0", 20,0,1);
    TH1F* h_taub1 = new TH1F("h_taub1", "h_taub1", 20,0,1);
    TH1F* h_taua1 = new TH1F("h_taua1", "h_taua1", 20,0,1);
    
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
        //2.0 tight PF jet selection
        int nFATJets = data.GetInt("FATnJet");
        if (nFATJets < 2) continue;
        vector<bool> &passFatJetTightID = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
        if (!passFatJetTightID[0]) continue;
        if (!passFatJetTightID[1]) continue;
        nPass[2]++;
        
        TClonesArray    *fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
        nPass[3]++;
        
        
        // take Leading and Trailing jet and observe HT
        float *jetMass = data.GetPtrFloat("FATjetPuppiSDmass");
        TLorentzVector* hzpJet = new TLorentzVector();
        Float_t HT = 0;
        Float_t jetEta[2], jetPt[2], jetM[2];
        int hm = 0, hPt = 0;
        for(int ij=0; ij<nFATJets; ij++)
        {
            TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
            HT += thisJet->Pt();
            if (ij>1) continue;
            if (ij==0) h_ptb->Fill(thisJet->Pt());
            if (ij==1) h_ptb1->Fill(thisJet->Pt());
            if (fabs(thisJet->Eta())>2.4) continue;
            if (thisJet->Pt()>300) hPt++;
            if (jetMass[ij]>105 && jetMass[ij]<135) hm++;
   
        }
        if (hPt<1) continue;
        nPass[5]++;
        TLorentzVector* thisjet = (TLorentzVector*)fatjetP4->At(0); 
        TLorentzVector* thatjet = (TLorentzVector*)fatjetP4->At(1);
        h_pta->Fill(thisjet->Pt());
        h_pta1->Fill(thatjet->Pt());
        
        // select fatjet 
        float   *FatjetTau1 = data.GetPtrFloat("FATjetPuppiTau1");
        float   *FatjetTau2 = data.GetPtrFloat("FATjetPuppiTau2");
        float Tau21[2];
        for (int i=0;i<2;i++) {
            Tau21[i]= FatjetTau2[i] / FatjetTau1[i];
        }
        h_taub->Fill(Tau21[0]);
        h_taub1->Fill(Tau21[1]);
        if (Tau21[0] > 0.55) continue;
        if (Tau21[1] > 0.55) continue;
        nPass[4]++;
        h_taua->Fill(Tau21[0]);
        h_taua1->Fill(Tau21[1]);
        
        
        if (hm<1) continue;
        nPass[6]++;
        if (fabs(jetEta[0]-jetEta[1]) > 1.3) continue; // Eta of LO and NLO fatjet
        nPass[7]++;
        //*hzpJet = *jetLO + *jetNLO;

        h_HT->Fill(HT);
        //h_4bJetM->Fill(hzpJet->M());
        for (int i=0;i<2;i++) {
            h_jetM[i]->Fill(jetM[i]);
            h_jetPt[i]->Fill(jetPt[i]);
            h_jetEta[i]->Fill(jetEta[i]);
        }
    } // end of loop over entries

    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    for(int i=0;i<20;i++) if(nPass[i]>0) std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;
    efferr(nPass[7],nTotal);
    
    // write efficiency
    string outputPath = "efficiency/" + outputFile;
    /*
    fstream file;
    file.open(outputPath,ios::out|ios::trunc);
    file << cutName[0] << "\t" << nTotal << "\n";  
    for (int i=0;i<7;i++) {
        file << cutName[i+1] << "\t" << float(nPass[i])/nTotal << "\n";
    }
    file.close();
    */
    string outputRootFile=Form("nTuple/%s_boost_w.root",zpmass.Data());
    //string outputRootFile=Form("%s_boost.root","1500");
    TFile* outFile = new TFile(outputRootFile.data(),"recreate");
    h_4bJetM->Write();
    h_HT->Write();
    for (int i=0;i<2;i++) {
        h_jetM[i]->Write();
        h_jetPt[i]->Write();
        h_jetEta[i]->Write();
        h_taub->Write();
        h_taub1->Write();
        h_taua->Write();
        h_taua1->Write();

    }
    outFile->Close();
    
}
