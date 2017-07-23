
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
void boosted_genmatch(int w, std::string inputFile){
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
        //2.0 tight PF jet selection
        int nFATJets = data.GetInt("FATnJet");
        if (nFATJets < 2) continue;
        vector<bool> &passFatJetTightID = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
        if (!passFatJetTightID[0]) continue;
        if (!passFatJetTightID[1]) continue;
        nPass[2]++;
        
        TClonesArray    *fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
        TClonesArray    *genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
        int    *genMomParId = data.GetPtrInt("genMomParId");
        nPass[3]++;
        
        // take genPar and find the higgs and zp jet, pdg=25 or 9000001
        int *genParId = data.GetPtrInt("genParId");
        int *da1 = data.GetPtrInt("genDa1");
        int *da2 = data.GetPtrInt("genDa2");
        int Hindex = -1;
        int zpindex[2] = {-1,-1};
        int HdaIndex[2] = {-1,-1};
        for (int i=0 ; i<30 ; i++) {
            if (genMomParId[i]==9000001&&zpindex[0]<0) zpindex[0] = i;
            else if (genMomParId[i]==9000001&&zpindex[1]<0) zpindex[1] = i;
            if (genParId[i]==25) {
                Hindex = i;
            }
            if (genMomParId[i]==25&&HdaIndex[0]<0) HdaIndex[0] = i;
            else if (genMomParId[i]==25&&HdaIndex[1]<0) HdaIndex[1] = i;

            //if (Hindex>=0&&zpindex[1]>=0) break;
        }
        if (Hindex<0||zpindex[1]<0||HdaIndex[1]<0) continue;
        nPass[4]++;
        TLorentzVector* genHJet = (TLorentzVector*)genParP4->At(Hindex);
        TLorentzVector* genZptobJet0 = (TLorentzVector*)genParP4->At(zpindex[0]);
        TLorentzVector* genZptobJet1 = (TLorentzVector*)genParP4->At(zpindex[1]);
        TLorentzVector* genHtobJet0 = (TLorentzVector*)genParP4->At(HdaIndex[0]);
        TLorentzVector* genHtobJet1 = (TLorentzVector*)genParP4->At(HdaIndex[1]);
        TLorentzVector* ZptobbJet = new TLorentzVector();
        *ZptobbJet = *genZptobJet0 + *genZptobJet1;
        h_zpGenDeltaR->Fill(genZptobJet0->DeltaR(*genZptobJet1));
        h_hGenDeltaR->Fill(genHtobJet0->DeltaR(*genHtobJet1));
        // take Leading and Trailing jet and observe HT
        float *jetSDmass = data.GetPtrFloat("FATjetPuppiSDmass");
        float   *FatjetTau1 = data.GetPtrFloat("FATjetPuppiTau1");
        float   *FatjetTau2 = data.GetPtrFloat("FATjetPuppiTau2");
        float Tau21[2];
        int matchHIndex = -1, matchzpIndex = -1;
        
        // Higgs
        for(int ij=0; ij<nFATJets; ij++)
        {
            
            TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
            if (thisJet->DeltaR(*genHJet)<0.8) {
                matchHIndex=ij;
                break;
            }
        }
        //if (matchHIndex<0 ) continue;
        nPass[5]++;
        
        if (matchHIndex>=0) {
            TLorentzVector* hjet = (TLorentzVector*)fatjetP4->At(matchHIndex); 
            Tau21[0]= FatjetTau2[matchHIndex] / FatjetTau1[matchHIndex];
            h_hM->Fill(jetSDmass[matchHIndex]);
            h_hpt->Fill(hjet->Pt());
            h_heta->Fill(hjet->Eta());
            h_htau->Fill(Tau21[0]);
        }
        
        // Zp
        for(int ij=0; ij<nFATJets; ij++) {
            TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
            //if (thisJet->DeltaR(*genZptobJet0)<0.8 && thisJet->DeltaR(*genZptobJet1)<0.8) {
            if (thisJet->DeltaR(*ZptobbJet)<0.8) {
                matchzpIndex = ij;
                break;
            }
        }
        //if (matchzpIndex<0 || matchHIndex==matchzpIndex) continue;
        
        nPass[6]++;
        if (matchzpIndex>=0) {
            TLorentzVector* zpjet = (TLorentzVector*)fatjetP4->At(matchzpIndex);
            Tau21[1]= FatjetTau2[matchzpIndex] / FatjetTau1[matchzpIndex];
            h_zpM->Fill(jetSDmass[matchzpIndex]);
            h_zppt->Fill(zpjet->Pt());
            h_zpeta->Fill(zpjet->Eta());
            h_zptau->Fill(Tau21[1]);
        }
        
    } // end of loop over entries

    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    for(int i=0;i<20;i++) if(nPass[i]>0) std::cout << "nPass[" << i << "] = " << nPass[i] << std::endl;
    efferr(nPass[6],nTotal);
    
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
    
}
