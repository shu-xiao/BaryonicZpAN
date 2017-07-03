
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
void boosted_xAna_BZ(std::string inputFile){

    //get TTree from file ...
    TreeReader data(inputFile.data());

    TString zpmass;
    zpmass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*_bb_MZp}; test=${test1%%_gq*}; echo \"${test}\"",inputFile.data()));
    string outputFile=Form("effiMZp%s.txt",zpmass.Data());

    TCanvas* c1 = new TCanvas("c1","",889*1.5,768);
    TH1F* h_ZpPt = new TH1F("h_ZpPt", "ZpPt", 25,0,2500);
    TH1F* h_ZpEta = new TH1F("h_ZpEta", "ZpEta", 30,-3,3);
    TH1F* h_ZpJetM = new TH1F("h_ZpJetM", "ZpJetM", 30,0,4500);
    TH1F* h_higgsPt = new TH1F("h_higgsPt", "higgs Pt", 20,0,2000);
    TH1F* h_higgsEta = new TH1F("h_higgsEta", "higgs Eta", 30,-3,3);
    TH1F* h_higgsJetM = new TH1F("h_higgsJetMass", "higgs jet Mass", 25,100,150);
    TH1F* h_4bJetM = new TH1F("h_4bMass", "4b mass", 16,100,3300);
    TH1F* h_HT = new TH1F("h_HT", "HT", 20,0,1500);

    Int_t nPass[20]={0};
    int Hindex[2]={-1,-1};
    string cutName[6] = {"total","goodVtx","\\\\trig","PassID", "findJet", "diJetCut"};
    //for(Long64_t jEntry=0; jEntry<100 ;jEntry++){
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

        if (jEntry % 20000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

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
            if ( tri[0] && tri[1] && tri[2] && tri[3] && tri[4] && tri[5] && tri[6] && tri[7] && results==1 )    
                {
                    //cout << thisTrig << endl;
                    passTrigger=true;
                    break;
                }
        } // end of trigger

        //if(!passTrigger) continue;
        nPass[1]++;
        //2.0 tight PF jet selection
        int nFATJets = data.GetInt("FATnJet");
        if (nFATJets < 2) continue;
        vector<bool> &passFatJetTightID = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
        //if (!passFatJetTightID[0]) continue;

        Float_t* neuHadFrac = data.GetPtrFloat("FATjetNHadEF");
        Float_t* neuEmFrac = data.GetPtrFloat("FATjetNEmEF");
        //Float_t* neuCons = data.GetPtrFloat("FATjet_nSV");
        Float_t* muoFrac = data.GetPtrFloat("FATjetMuoEF");
        Float_t* chaHadFrac = data.GetPtrFloat("FATjetCHadEF");
        Int_t* chaMulti = data.GetPtrInt("FATjetCMulti");
        Float_t* chaEmFrac = data.GetPtrFloat("FATjetCEmEF");
        
        if (neuHadFrac[0] > 0.9) continue;
        if (neuEmFrac[0] > 0.9)  continue;
        ////if (neuCons[0] <= 1) continue;
        if (muoFrac[0] > 0.8) continue;
        if (chaHadFrac[0] <= 0) continue;
        if (chaMulti[0] <= 0) continue;
        if (chaEmFrac[0] > 0.9) continue;
        
        // select fatjet 

        TClonesArray    *fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
        TClonesArray    *genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
        float   *FatjetTau1 = data.GetPtrFloat("FATjetPuppiTau1");
        float   *FatjetTau2 = data.GetPtrFloat("FATjetPuppiTau2");
        float Tau21 = *FatjetTau2 / *FatjetTau1;
     
        vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("FATsubjetSDPx", nFATJets);
        vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("FATsubjetSDPy", nFATJets); // integral number need fix
        vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("FATsubjetSDPz", nFATJets);
        vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("FATsubjetSDE", nFATJets);  // subjet 4 vector
        
        // take LO and NLO jet and observe HT
        /*
        int HIndex[2]={-1,-1};
        Float_t HT = 0;
        Float_t jetEta[2];
        for(int ij=0; ij<nFATJets; ij++)
        {
            TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
            if (ij<2) jetEta[ij] = thisJet->Eta();
            HT += thisJet->Pt();
            if(thisJet->Pt()<300)continue;
            if(fabs(thisJet->Eta())>2.4)continue;
            if(!passFatJetTightID[ij])continue;
            //if (HIndex[0]<0) HIndex[0]=ij;
            //else if (Hindex[1]<0) HIndex[1]=ij;
            //else break;
        }
        h_HT->Fill(HT);
        //if(HIndex[1]<0) continue;
        if (fabs(jetEta[0]-jetEta[1]) > 1.3) continue; // Eta of LO and NLO fatjet
        */
        nPass[2]++;
        
        // Jet kinetic selection
        // select b quarks
        TLorentzVector* subjetP4[2][2];
        for(int i=0;i<2;i++){
            //if( !passFatJetTightID[i] )continue;
            for(int j=0;j<2;j++){
                subjetP4[i][j]=new TLorentzVector();
                //subjetP4[i][j]->SetPxPyPzE(subjetSDPx[HIndex[i]][j],subjetSDPy[HIndex[i]][j],subjetSDPz[HIndex[i]][j],subjetSDE[HIndex[i]][j]);
            }
        }
        nPass[3]++; 
        // Scan for DM, 9000001
        Int_t *parID = data.GetPtrInt("genParId");
        Int_t *genDa1 = data.GetPtrInt("genDa1");
        Int_t *genDa2 = data.GetPtrInt("genDa2");
        Int_t *genMomParId = data.GetPtrInt("genMomParId");
        Int_t *genParSt = data.GetPtrInt("genParSt");
        Int_t nPar = data.GetInt("nGenPar");
        Int_t zpIndex[2] = {-1}, higgsIndex[2] = {-1};
        Int_t nZpmo = 0, nHiggsmo = 0;
        for (int ij=0;ij<30;ij++) {
            if (genParSt[ij] != 23) continue;
            if (genMomParId[ij]==9000001) {
                nZpmo += 1;
                if (zpIndex[0] < 0) zpIndex[0] = ij;
                else if (zpIndex[1] < 0) zpIndex[1] = ij;
            }
            if (genMomParId[ij]==25) {
                nHiggsmo += 1;
                if (higgsIndex[0] < 0) higgsIndex[0] = ij;
                else if (higgsIndex[1] < 1) higgsIndex[1] = ij;
            }
            if (higgsIndex[1] > 0 && zpIndex[1] > 0) break;
        }
        if (nZpmo < 2 || nHiggsmo < 2) continue;
        
        TLorentzVector* zpJetDa1 = (TLorentzVector*)genParP4->At(zpIndex[0]);
        TLorentzVector* zpJetDa2 = (TLorentzVector*)genParP4->At(zpIndex[1]);
        TLorentzVector* zpJet = new TLorentzVector();
        *zpJet = *zpJetDa1 + *zpJetDa2;
        TLorentzVector* higgsJetDa1 = (TLorentzVector*)genParP4->At(higgsIndex[0]);
        TLorentzVector* higgsJetDa2 = (TLorentzVector*)genParP4->At(higgsIndex[1]);
        TLorentzVector* higgsJet = new TLorentzVector();
        *higgsJet = *higgsJetDa1 + *higgsJetDa2;
        TLorentzVector* b4Jet = new TLorentzVector();
        *b4Jet = *higgsJet + *zpJet;
        float b4Mass = b4Jet->M();
        float zpJetM = zpJet->M();
        float zpJetPt = zpJet->Pt();
        float zpJetEta = zpJet->Eta();
        float higgsJetM = higgsJet->M();
        float higgsJetPt = higgsJet->Pt();
        float higgsJetEta = higgsJet->Eta();
        
        Float_t*  fatjetPuppiSDmass = data.GetPtrFloat("FATjetPuppiSDmass");
        //if ((fatjetPuppiSDmass[0] < 100 || fatjetPuppiSDmass[0] > 150) && (fatjetPuppiSDmass[1] < 100 || fatjetPuppiSDmass[1] > 150)) continue;
        nPass[4]++;
        // assume MZp > M_higss
        /*
        TLorentzVector* higgsJet = (TLorentzVector*)fatjetP4->At(HIndex[0]);
        TLorentzVector* ZpJet = (TLorentzVector*)fatjetP4->At(HIndex[1]);
        TLorentzVector *diJet = new TLorentzVector();
        //cout << "HIndex=0: " << higgsJet->M() << "\tHIndex=1: " << ZpJet->M() << endl;
        
        if (higgsJet->M() > ZpJet->M())
        {
            higgsJet = (TLorentzVector*)fatjetP4->At(HIndex[1]);
            ZpJet = (TLorentzVector*)fatjetP4->At(HIndex[0]);
        }
        
        //crab_50jobs/crab_MonoH-_ZpBaryonic_bb_MZp900_gq25_MDM100_hbbcout << "Index1: " << fatjetPuppiSDmass[Hindex[0]] << "\tIndex2: " <<fatjetPuppiSDmass[Hindex[1]] << endl;
        // mass loose cut

        // if (!((higgsJet->M() > 100 && higgsJet->M() < 150) || (ZpJet->M() > 100 && ZpJet->M() < 150))) continue;
        //if ((fatjetPuppiSDmass[Hindex[0]] < 100 || fatjetPuppiSDmass[Hindex[0]] > 150) && (fatjetPuppiSDmass[Hindex[1]] < 100 || fatjetPuppiSDmass[Hindex[1]] > 150)) continue;
        *diJet = *higgsJet + *ZpJet;
       
        float redMass2Jet = diJet->M() - (higgsJet->M()-125) - (ZpJet->M()-std::stof(zpmass.Data()));
        //if (redMass2Jet < 200) continue;
        float higgsPt, higgsEta, ZpPt, ZpEta;
        higgsPt = higgsJet->Pt();
        higgsEta = higgsJet->Eta();
        ZpPt = ZpJet->Pt();
        ZpEta = ZpJet->Eta();
        // if (abs(higgsEta-ZpEta) > 1.3) continue;
        */
        float redMass2Jet = b4Jet->M() - ( higgsJetM -125 ) - ( zpJetM - 1500 );
        if (redMass2Jet < 200) continue;
        if (abs(zpJetEta-higgsJetEta) > 1.3) continue;
        if (Tau21 > 0.55) continue;
        nPass[5]++;
        
        h_higgsJetM->Fill(higgsJetM);
        h_higgsPt->Fill(higgsJetPt);
        h_higgsEta->Fill(higgsJetEta);
        h_ZpJetM->Fill(zpJetM);    
        h_ZpPt->Fill(zpJetPt);
        h_ZpEta->Fill(zpJetEta);
        h_4bJetM->Fill(b4Mass);
        //h_diJetM->Fill(redMass2Jet); 
    
    } // end of loop over entries

    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    for(int i=0;i<20;i++) if(nPass[i]>0) std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;
    efferr(nPass[5],nTotal);
    
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
    h_higgsPt->Draw("hist");
    c1->SaveAs("higgsPt.png");
    h_higgsEta->Draw("hist");
    c1->SaveAs("higgsEta.png");
    h_higgsJetM->Draw("hist");
    c1->SaveAs("higgsJetM.png");
    h_ZpPt->Draw("hist");
    c1->SaveAs("ZpPt.png");
    h_ZpEta->Draw("hist");
    c1->SaveAs("ZpEta.png");
    h_ZpJetM->Draw("hist");
    c1->SaveAs("ZpJetM.png");
    h_4bJetM->Draw("hist");
    c1->SaveAs("4bMass.png");
    //h_HT->Draw("hist");
    //c1->SaveAs("HT.png");
    
    /*
    outputRootFile=Form("%s_boost.root",zpmass.Data());
    TFile* outFile = new TFile(outputRootFile,"recreate");
    outFile->Close();
    */
}

