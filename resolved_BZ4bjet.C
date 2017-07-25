
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

    TH1F* h_htobbpt0 = new TH1F("h_htobbpt0","h_htobbpt0_match",20,0,1400);
    TH1F* h_htobbpt1 = new TH1F("h_htobbpt1","h_htobbpt1_match",20,0,1400);
    TH1F* h_zptobbpt0 = new TH1F("h_zptobbpt0","h_zptobbpt0_match",20,0,1400);
    TH1F* h_zptobbpt1 = new TH1F("h_zptobbpt1","h_zptobbpt1_match",20,0,1400);

    TH1F* h_hCIS0 = new TH1F("h_higgsCIS0", "h_higgstobCISVV2_0", 20,0,1);
    TH1F* h_hCIS1 = new TH1F("h_higgsCIS1", "h_higgstobCISVV2_1", 20,0,1);
    TH1F* h_zpCIS0 = new TH1F("h_zpCIS0", "h_zptobCISVV2_0", 20,0,1);
    TH1F* h_zpCIS1 = new TH1F("h_zpCIS1", "h_zptobCISVV2_1", 20,0,1);
    
    TH1F* h_hM = new TH1F("h_higgsM", "h_higgsM_match", 25,0,600);
    TH1F* h_zpM = new TH1F("h_ZpM", "h_ZpM_match", 25,0,1800);
    
    TH1F* h_hDeltaR = new TH1F("h_HiggstobbdeltaR", "h_HiggstobbDeltaR_match", 25,0,6);
    TH1F* h_zpDeltaR = new TH1F("h_ZptobbDeltaR", "h_ZptobbDeltaR_match", 25,0,6);
    
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
        float* thinJetCSV =  data.GetPtrFloat("THINjetCISVV2");
        TClonesArray    *genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
        int    *genMomParId = data.GetPtrInt("genMomParId");
        
        // take genPar and find the higgs and zp jet, pdg=25 or 9000001
        int *genParId = data.GetPtrInt("genParId");
        int *da1 = data.GetPtrInt("genDa1");
        int *da2 = data.GetPtrInt("genDa2");
        int Hgenindex = -1;
        int zpdaGenIndex[2] = {-1,-1};
        int HdaGenIndex[2] = {-1,-1};
        for (int i=0 ; i<30 ; i++) {
            if (genMomParId[i]==9000001&&zpdaGenIndex[0]<0) zpdaGenIndex[0] = i;
            else if (genMomParId[i]==9000001&&zpdaGenIndex[1]<0) zpdaGenIndex[1] = i;
            if (genParId[i]==25) {
                Hgenindex = i;
            }
            if (genMomParId[i]==25&&HdaGenIndex[0]<0) HdaGenIndex[0] = i;
            else if (genMomParId[i]==25&&HdaGenIndex[1]<0) HdaGenIndex[1] = i;

            //if (Hgenindex>=0&&zpdaGenIndex[1]>=0) break;
        }
        if (Hgenindex<0||zpdaGenIndex[1]<0||HdaGenIndex[1]<0) continue;
        nPass[2]++;
        
        // genMatch
        TLorentzVector* genHJet = (TLorentzVector*)genParP4->At(Hgenindex);
        TLorentzVector* genZptobJet0 = (TLorentzVector*)genParP4->At(zpdaGenIndex[0]);
        TLorentzVector* genZptobJet1 = (TLorentzVector*)genParP4->At(zpdaGenIndex[1]);
        TLorentzVector* genHtobJet0 = (TLorentzVector*)genParP4->At(HdaGenIndex[0]);
        TLorentzVector* genHtobJet1 = (TLorentzVector*)genParP4->At(HdaGenIndex[1]);
        const int nTHINJets     = data.GetInt("THINnJet");
        TClonesArray* thinjetP4 = (TClonesArray*) data.GetPtrTObject("THINjetP4");
        int matchHIndex[2] = {-1,-1}, matchZpIndex[2] = {-1,-1};

        for(int k=0; k<nTHINJets; k++)
        {
            TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(k);
            if (thisJet->DeltaR(*genHtobJet0)<0.4 && matchHIndex[0]<0) {
                matchHIndex[0]=k;
                continue;
            }
            if (thisJet->DeltaR(*genHtobJet1)<0.4 && matchHIndex[1]<0) {
                matchHIndex[1]=k;
                continue;
            }
            if (thisJet->DeltaR(*genZptobJet0)<0.4 && matchZpIndex[0]<0) {
                matchZpIndex[0]=k;
                continue;
            }
            if (thisJet->DeltaR(*genZptobJet1)<0.4 && matchZpIndex[1]<0) {
                matchZpIndex[1]=k;
                continue;
            }
        }
        if (matchHIndex[0]<0||matchHIndex[1]<0||matchZpIndex[0]<0||matchZpIndex[1]<0) continue;
        TLorentzVector* matchHtobJet0 = (TLorentzVector*)thinjetP4->At(matchHIndex[0]);
        TLorentzVector* matchHtobJet1 = (TLorentzVector*)thinjetP4->At(matchHIndex[1]);
        TLorentzVector* matchZptobJet0 = (TLorentzVector*)thinjetP4->At(matchZpIndex[0]);
        TLorentzVector* matchZptobJet1 = (TLorentzVector*)thinjetP4->At(matchZpIndex[1]);

        h_zpM->Fill((*matchZptobJet0+*matchZptobJet1).M()); 
        h_zppt->Fill((*matchZptobJet0+*matchZptobJet1).Pt());
        h_hM->Fill((*matchHtobJet0+*matchHtobJet1).M());
        h_hpt->Fill((*matchHtobJet0+*matchHtobJet1).Pt());
        h_htobbpt0->Fill(matchHtobJet0->Pt());
        h_htobbpt1->Fill(matchHtobJet1->Pt());
        h_zptobbpt0->Fill(matchZptobJet0->Pt());
        h_zptobbpt1->Fill(matchZptobJet1->Pt());

        h_hCIS0->Fill(thinJetCSV[matchHIndex[0]]);
        h_hCIS1->Fill(thinJetCSV[matchHIndex[1]]);
        h_zpCIS0->Fill(thinJetCSV[matchZpIndex[0]]);
        h_zpCIS1->Fill(thinJetCSV[matchZpIndex[1]]);

        h_hDeltaR->Fill(matchHtobJet0->DeltaR(*matchHtobJet1));
        h_zpDeltaR->Fill(matchZptobJet0->DeltaR(*matchZptobJet1));

        // reco higgs and Z'
        if (nTHINJets<4) continue;
        vector<bool>& passThinJetLooseID = *((vector<bool>*) data.GetPtr("THINjetPassIDLoose"));
        vector<bool>& passThinJetPUID = *((vector<bool>*) data.GetPtr("THINisPUJetIDMedium"));    
    
        int Hindex[2]={-1,-1};
        int Zpindex[2]={-1,-1};
        float maxHpt=-999;
        float maxZppt=-999;
        const float ptCut = 30;
        const float thinCISVV2cut = 0;
        //const float thinCISVV2cut = 0.5426;
        for(int ij=0; ij < nTHINJets; ij++){
            TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(ij);
            if(thisJet->Pt()<ptCut)continue;
            if(fabs(thisJet->Eta())>2.4)continue;
            if(!passThinJetLooseID[ij])continue;
            if(!passThinJetPUID[ij])continue;
      
        // for b-jet (medium ID)
            if(thinJetCSV[ij]<thinCISVV2cut)continue;

            for(int jj=0; jj < ij ; jj++){
	        TLorentzVector* thatJet = (TLorentzVector*)thinjetP4->At(jj);
	        if(thatJet->Pt()<ptCut)continue;
	        if(fabs(thatJet->Eta())>2.4)continue;
	        if(!passThinJetLooseID[jj])continue;
	        if(!passThinJetPUID[jj])continue;
	
	// for b-jet (medium ID)
	        if(thinJetCSV[jj]<thinCISVV2cut)continue;

	        float thisHpt = (*thisJet + *thatJet).Pt();
	        float thisHmass = (*thisJet + *thatJet).M();
	        //if(thisHpt < 150.)continue;

	        //if(thisHmass < 100)continue;
	        //if(thisHmass > 150)continue;
	        if(thisHpt>maxHpt)
	        {
	            Hindex[0]=jj;
	            Hindex[1]=ij;
	            maxHpt   =thisHpt;
	        }

            } // end of inner loop jet


        } // end of outer loop jet


        if(Hindex[0]<0 || Hindex[1]<0)continue;
        nPass[3]++;



        // Zp
        for(int ij=0; ij < nTHINJets; ij++){
            if (ij==Hindex[1]) continue;
            TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(ij);
            if(thisJet->Pt()<ptCut)continue;
            if(fabs(thisJet->Eta())>2.4)continue;
            if(!passThinJetLooseID[ij])continue;
            if(!passThinJetPUID[ij])continue;
      
      // for Zp-jet (medium ID)
            if(thinJetCSV[ij]<thinCISVV2cut)continue;

            for(int jj=0; jj < ij ; jj++){
                if (jj==Hindex[0]) continue;
                TLorentzVector* thatJet = (TLorentzVector*)thinjetP4->At(jj);
	        if(thatJet->Pt()<ptCut)continue;
	        if(fabs(thatJet->Eta())>2.4)continue;
	        if(!passThinJetLooseID[jj])continue;
	        if(!passThinJetPUID[jj])continue;
	
	// for Zp-jet (medium ID)
	        if(thinJetCSV[jj]<thinCISVV2cut)continue;

	        float thisZppt = (*thisJet + *thatJet).Pt();
	        float thisZpmass = (*thisJet + *thatJet).M();
	        //if(thisZppt < 150.)continue;

	        if(thisZppt>maxZppt)
	        {
	            Zpindex[0]=jj;
	            Zpindex[1]=ij;
	            maxZppt   =thisZppt;
	        }

            } // end of inner loop jet


        } // end of outer loop jet


    if(Zpindex[0]<0 || Zpindex[1]<0)continue;
    nPass[4]++;
    } // end of loop over entries
    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    for(int i=0;i<20;i++) if(nPass[i]>0) std::cout << "nPass[" << i << "] = " << nPass[i] << std::endl;
    efferr(nPass[4],nTotal);
    
    
    h_zpM->Draw("hist");
    c1->Print("thinJetgenMatch.pdf(");
    h_hM->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_zppt->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_hpt->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_htobbpt0->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_htobbpt1->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_zptobbpt0->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_zptobbpt1->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_hCIS0->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_hCIS1->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_zpCIS0->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_zpCIS1->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_hDeltaR->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_zpDeltaR->Draw("hist");
    c1->Print("thinJetgenMatch.pdf)");




/*
    // write efficiency
    string outputPath = "efficiency/" + outputFile;
    string outputRootFile=Form("nTuple/%s_boost_match.root",zpmass.Data());
    
    TFile* outFile = new TFile(outputRootFile.data(),"recreate");
    h_hM->Write();
    h_zpM->Write();
    h_hpt->Write();
    h_zppt->Write();
    outFile->Close();
   */ 
}
