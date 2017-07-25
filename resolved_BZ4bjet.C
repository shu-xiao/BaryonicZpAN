
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
void resolved_BZ4bjet(int w, std::string inputFile){
    setNCUStyle(true);
    //get TTree from file ...
    TreeReader data(inputFile.data());

    TString zpmass;
    zpmass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*_bb_MZp}; test=${test1%%.root}; echo \"${test}\"",inputFile.data()));
    string outputFile=Form("effiMZp%s.txt",zpmass.Data());

    TCanvas* c1 = new TCanvas("c1","",600,600);

    TH1F* h_hpt = new TH1F("h_higgsPt", "h_higgsPt_match", 20,0,1400);
    TH1F* h_zppt = new TH1F("h_ZpPt", "h_ZpPt_match", 20,0,1400);

    TH1F* h_htau = new TH1F("h_higgsTau21", "h_higgsTau21_match", 20,0,1);
    TH1F* h_zptau = new TH1F("h_ZpTau21", "h_ZpTau21_match", 20,0,1);
    
    TH1F* h_heta = new TH1F("h_higgsEta", "h_higgsEta_match", 20,-3,3);
    TH1F* h_zpeta = new TH1F("h_ZpEta", "h_ZpEta_match", 20,-3,3);
    TH1F* h_hM = new TH1F("h_higgsM", "h_higgsM_match", 25,0,600);
    TH1F* h_zpM = new TH1F("h_ZpM", "h_ZpM_match", 25,0,600);
    Int_t nPass[20]={0};
    TH1F* h_hGenDeltaR = new TH1F("h_HiggstobbdeltaR", "h_HiggstobbDeltaR_gen", 25,0,6);
    TH1F* h_zpGenDeltaR = new TH1F("h_ZptobbDeltaR", "h_ZptobbDeltaR_gen", 25,0,6);
    

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
        
        const int nTHINJets     = data.GetInt("THINnJet");
        TClonesArray* thinjetP4 = (TClonesArray*) data.GetPtrTObject("THINjetP4");
        float* thinJetCSV =  data.GetPtrFloat("THINjetCISVV2");
        vector<bool>& passThinJetLooseID = *((vector<bool>*) data.GetPtr("THINjetPassIDLoose"));
        vector<bool>& passThinJetPUID = *((vector<bool>*) data.GetPtr("THINisPUJetIDMedium"));    
    
        int Hindex[2]={-1,-1};
        int Zpindex[2]={-1,-1};
        float maxHpt=-999;
        float maxZppt=-999;
        float ptCut = 30;
        float thinCISVV2 = 0.5426;
        for(int ij=0; ij < nTHINJets; ij++){
            TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(ij);
            if(thisJet->Pt()<ptCut)continue;
            if(fabs(thisJet->Eta())>2.4)continue;
            if(!passThinJetLooseID[ij])continue;
            if(!passThinJetPUID[ij])continue;
      
      // for b-jet (medium ID)
            if(thinJetCSV[ij]<thinCISVV2)continue;

            for(int jj=0; jj < ij ; jj++){
	        TLorentzVector* thatJet = (TLorentzVector*)thinjetP4->At(jj);
	        if(thatJet->Pt()<ptCut)continue;
	        if(fabs(thatJet->Eta())>2.4)continue;
	        if(!passThinJetLooseID[jj])continue;
	        if(!passThinJetPUID[jj])continue;
	
	// for b-jet (medium ID)
	        if(thinJetCSV[jj]<thinCISVV2)continue;

	        float thisHpt = (*thisJet + *thatJet).Pt();
	        float thisHmass = (*thisJet + *thatJet).M();
	        if(thisHpt < 150.)continue;

	        if(thisHmass < 100)continue;
	        if(thisHmass > 150)continue;
	
	        if(thisHpt>maxHpt)
	        {
	            Hindex[0]=jj;
	            Hindex[1]=ij;
	            maxHpt   =thisHpt;
	        }

            } // end of inner loop jet


        } // end of outer loop jet


        if(Hindex[0]<0 || Hindex[1]<0)continue;
        nPass[2]++;

        for(int ij=0; ij < nTHINJets; ij++){
            if (ij==Hindex[1]) continue;
            TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(ij);
            if(thisJet->Pt()<ptCut)continue;
            if(fabs(thisJet->Eta())>2.4)continue;
            if(!passThinJetLooseID[ij])continue;
            if(!passThinJetPUID[ij])continue;
      
      // for b-jet (medium ID)
            if(thinJetCSV[ij]<thinCISVV2)continue;

            for(int jj=0; jj < ij ; jj++){
                if (jj==Hindex[0]) continue;
                TLorentzVector* thatJet = (TLorentzVector*)thinjetP4->At(jj);
	        if(thatJet->Pt()<ptCut)continue;
	        if(fabs(thatJet->Eta())>2.4)continue;
	        if(!passThinJetLooseID[jj])continue;
	        if(!passThinJetPUID[jj])continue;
	
	// for b-jet (medium ID)
	        if(thinJetCSV[jj]<thinCISVV2)continue;

	        float thisZppt = (*thisJet + *thatJet).Pt();
	        float thisZpmass = (*thisJet + *thatJet).M();
	        //if(thisZppt < 150.)continue;

	
	        if(thisZppt>maxZppt)
	        {
	            Hindex[0]=jj;
	            Hindex[1]=ij;
	            maxZppt   =thisZppt;
	        }

            } // end of inner loop jet


        } // end of outer loop jet


    if(Zpindex[0]<0 || Zpindex[1]<0)continue;
    nPass[3]++;
    } // end of loop over entries
    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    for(int i=0;i<20;i++) if(nPass[i]>0) std::cout << "nPass[" << i << "] = " << nPass[i] << std::endl;
    efferr(nPass[3],nTotal);
    
    /*
    h_zpM->Draw("hist");
    c1->Print("genMatch.pdf(");
    h_hM->Draw("hist");
    c1->Print("genMatch.pdf");
    h_zppt->Draw("hist");
    c1->Print("genMatch.pdf");
    h_hpt->Draw("hist");
    c1->Print("genMatch.pdf");
    h_zptau->Draw("hist");
    c1->Print("genMatch.pdf");
    h_htau->Draw("hist");
    c1->Print("genMatch.pdf");
    h_heta->Draw("hist");
    c1->Print("genMatch.pdf");
    h_zpeta->Draw("hist");
    c1->Print("genMatch.pdf");
    h_hGenDeltaR->Draw("hist");
    c1->Print("genMatch.pdf");
    h_zpGenDeltaR->Draw("hist");
    c1->Print("genMatch.pdf)");





    // write efficiency
    string outputPath = "efficiency/" + outputFile;
    string outputRootFile=Form("nTuple/%s_boost_match.root",zpmass.Data());
    
    TFile* outFile = new TFile(outputRootFile.data(),"recreate");
    h_hM->Write();
    h_zpM->Write();
    h_hpt->Write();
    h_zppt->Write();
    h_htau->Write();
    h_zptau->Write();
    h_heta->Write();
    h_zpeta->Write();
    h_zpGenDeltaR->Write();
    h_hGenDeltaR->Write();
    outFile->Close();
   */ 
}
