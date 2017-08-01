
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
#include "TMath.h"
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
void boosted_bb2HDMgenmatch(int w, std::string inputFile){
    setNCUStyle(true);
    //get TTree from file ...
    TreeReader data(inputFile.data());

    TString a0mass;
    a0mass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*_bb_MA0}; test=${test1%%.root}; echo \"${test}\"",inputFile.data()));
    string outputFile=Form("effiMA0%s.txt",a0mass.Data());

    TCanvas* c1 = new TCanvas("c1","",600,600);

    TH1F* h_hpt = new TH1F("h_higgsPt", "h_higgsPt_match", 20,0,1400);
    TH1F* h_a0pt = new TH1F("h_A0Pt", "h_A0Pt_match", 20,0,1400);

    TH1F* h_htau = new TH1F("h_higgsTau21", "h_higgsTau21_match", 20,0,1);
    TH1F* h_a0tau = new TH1F("h_A0Tau21", "h_A0Tau21_match", 20,0,1);
    
    TH1F* h_heta = new TH1F("h_higgsEta", "h_higgsEta_match", 20,-3,3);
    TH1F* h_a0eta = new TH1F("h_A0Eta", "h_A0Eta_match", 20,-3,3);
    
    TH1F* h_hM = new TH1F("h_higgsM", "h_higgsM_match", 25,0,600);
    TH1F* h_a0M = new TH1F("h_A0M", "h_A0M_match", 25,0,600);
    
    TH1F* h_nCA15jet = new TH1F("h_nCA15jet", "h_nCA15jet", 6,-0.5,5.5);
    
    TH1F* h_hbbtag = new TH1F("h_higgsdoublebtag","h_higgsdoublebtag_match",40,-1,1);
    TH1F* h_a0bbtag = new TH1F("h_a0doublebtag","h_a0doublebtag_match",40,-1,1);

    TH1F* h_hCISVV2 = new TH1F("h_higgsCISVV2","h_higgsCISVV2_match",20,0,1);
    TH1F* h_a0CISVV2 = new TH1F("h_a0CISVV2","h_a0CISVV2_match",20,0,1);

    TH1F* h_hECF_2_3_10 = new TH1F("h_higgsECF_2_3_10","h_higgsECF_2_3_10_match",50,0,0.1);
    TH1F* h_hECF_1_2_10 = new TH1F("h_higgsECF_1_2_10","h_higgsECF_1_2_10_match",50,0,0.5);
    TH1F* h_a0ECF_2_3_10 = new TH1F("h_a0ECF_2_3_10","h_a0ECF_2_3_10_match",50,0,0.1);
    TH1F* h_a0ECF_1_2_10 = new TH1F("h_a0ECF_1_2_10","h_a0ECF_1_2_10_match",50,0,0.5);
    TH1F* h_hN2 = new TH1F("h_hN2","h_hN2",50,0,1);
    TH1F* h_a0N2 = new TH1F("h_a0N2","h_a0N2",50,0,1);

    Int_t nPass[20]={0};
    //TH1F* h_hGenDeltaR = new TH1F("h_HiggstobbdeltaR", "h_HiggstobbDeltaR_gen", 25,0,6);
    //TH1F* h_a0GenDeltaR = new TH1F("h_A0tobbDeltaR", "h_A0tobbDeltaR_gen", 25,0,6);
    TH1F* h_a0HDeltaR_gen = new TH1F("h_A0HDeltaR_gen", "h_A0andhiggsDeltaR_gen", 25,0,6);
    TH1F* h_a0HDeltaR_match = new TH1F("h_A0HDeltaR_match", "h_A0andhiggsDeltaR_match", 25,0,6);
    

    //for(Long64_t jEntry=0; jEntry<1 ;jEntry++){
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        
        // broken events
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
        int nCA15Jets = data.GetInt("CA15PuppinJet");
        h_nCA15jet->Fill(nCA15Jets);
        if (nCA15Jets < 2) continue;
        vector<bool> &passCA15JetTightID = *((vector<bool>*) data.GetPtr("CA15PuppijetPassIDTight"));
        if (!passCA15JetTightID[0]) continue;
        if (!passCA15JetTightID[1]) continue;
        nPass[2]++;
        
        TClonesArray    *CA15jetP4 = (TClonesArray*) data.GetPtrTObject("CA15PuppijetP4");
        TClonesArray    *genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
        int    *genMomParId = data.GetPtrInt("genMomParId");
        nPass[3]++;
        
        // take genPar and find the higgs and a0 jet, pdg=25 or 28 
        int *genParId = data.GetPtrInt("genParId");
        int *da1 = data.GetPtrInt("genDa1");
        int *da2 = data.GetPtrInt("genDa2");
        int Hindex = -1;
        //int a0index[2] = {-1,-1};
        //int HdaIndex[2] = {-1,-1};
        int a0index = -1;
        for (int i=0 ; i<30 ; i++) {
            if (genParId[i]==28&&a0index<0) a0index = i;
            if (genParId[i]==25) {
                Hindex = i;
            }

            if (Hindex>=0&&a0index>=0) break;
        }
        if (Hindex<0||a0index<0) continue;
        nPass[4]++;
        TLorentzVector* genHJet = (TLorentzVector*)genParP4->At(Hindex);
        TLorentzVector* genA0Jet = (TLorentzVector*)genParP4->At(a0index);
        h_a0HDeltaR_gen->Fill(genHJet->DeltaR(*genA0Jet));
        // take Leading and Trailing jet and observe HT
        float *jetSDmass = data.GetPtrFloat("CA15PuppijetSDmass");
        float *CA15jetTau1 = data.GetPtrFloat("CA15PuppijetTau1");
        float *CA15jetTau2 = data.GetPtrFloat("CA15PuppijetTau2");
        float *CA15jetCISVV2 = data.GetPtrFloat("CA15PuppijetCISVV2");
        float *CA15jetDoublebtag = data.GetPtrFloat("CA15Puppi_doublebtag");
        float *CA15jetECF_2_3_10 = data.GetPtrFloat("CA15PuppiECF_2_3_10");
        float *CA15jetECF_1_2_10 = data.GetPtrFloat("CA15PuppiECF_1_2_10");
        float Tau21[2];
        int matchHIndex = -1, matcha0Index = -1;
        TLorentzVector *hjet, *a0jet;        
        // Higgs
        for(int ij=0; ij<nCA15Jets; ij++)
        {
            
            TLorentzVector* thisJet = (TLorentzVector*)CA15jetP4->At(ij);
            if (thisJet->DeltaR(*genHJet)<1.5) {
                matchHIndex=ij;
                break;
            }
        }
        if (matchHIndex<0 ) continue;
        nPass[5]++;
        
        if (matchHIndex>=0) {
            hjet = (TLorentzVector*)CA15jetP4->At(matchHIndex); 
            Tau21[0] = CA15jetTau2[matchHIndex] / CA15jetTau1[matchHIndex];
            h_hM->Fill(jetSDmass[matchHIndex]);
            h_hpt->Fill(hjet->Pt());
            h_heta->Fill(hjet->Eta());
            h_htau->Fill(Tau21[0]);
            h_hbbtag->Fill(CA15jetDoublebtag[matchHIndex]);
            h_hECF_2_3_10->Fill(CA15jetECF_2_3_10[matchHIndex]);
            h_hECF_1_2_10->Fill(CA15jetECF_1_2_10[matchHIndex]);
            h_hCISVV2->Fill(CA15jetCISVV2[matchHIndex]);
            float hN2 = CA15jetECF_2_3_10[matchHIndex]/pow(CA15jetECF_1_2_10[matchHIndex],2.00);
            h_hN2->Fill(hN2);
        }
        
        // A0
        for(int ij=0; ij<nCA15Jets; ij++) {
            if (ij==matchHIndex) continue;
            TLorentzVector* thisJet = (TLorentzVector*)CA15jetP4->At(ij);
            if (thisJet->DeltaR(*genA0Jet)<1.5) {
                matcha0Index = ij;
                break;
            }
        }
        if (matcha0Index<0) continue;
        
        nPass[6]++;
        if (matcha0Index>=0) {
            a0jet = (TLorentzVector*)CA15jetP4->At(matcha0Index);
            Tau21[1]= CA15jetTau2[matcha0Index] / CA15jetTau1[matcha0Index];
            h_a0M->Fill(jetSDmass[matcha0Index]);
            h_a0pt->Fill(a0jet->Pt());
            h_a0eta->Fill(a0jet->Eta());
            h_a0tau->Fill(Tau21[1]);
            h_a0bbtag->Fill(CA15jetDoublebtag[matcha0Index]);
            h_a0ECF_2_3_10->Fill(CA15jetECF_2_3_10[matcha0Index]);
            h_a0ECF_1_2_10->Fill(CA15jetECF_1_2_10[matcha0Index]);
            h_a0CISVV2->Fill(CA15jetCISVV2[matcha0Index]);
            float a0N2 = CA15jetECF_2_3_10[matcha0Index]/pow(CA15jetECF_1_2_10[matcha0Index],2.00);
            h_a0N2->Fill(a0N2);
        }
        if (matchHIndex>=0 && matcha0Index>=0) h_a0HDeltaR_match->Fill(a0jet->DeltaR(*hjet));
    } // end of loop over entries

    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    for(int i=0;i<20;i++) if(nPass[i]>0) std::cout << "nPass[" << i << "] = " << nPass[i] << std::endl;
    efferr(nPass[6],nTotal);
    
    h_nCA15jet->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf(");
    h_a0M->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_hM->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_a0pt->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_hpt->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_a0tau->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_htau->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_heta->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_a0eta->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_a0HDeltaR_gen->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_a0HDeltaR_match->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_hbbtag->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_a0bbtag->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_hECF_2_3_10->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_a0ECF_2_3_10->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_hECF_1_2_10->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_a0ECF_1_2_10->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_hN2->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_a0N2->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_hCISVV2->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf");
    h_a0CISVV2->Draw("hist");
    c1->Print("bb2HDM_MZp1500_MA0300.pdf)");





    // write efficiency
    string outputPath = "efficiency/" + outputFile;
    string outputRootFile="nTuple/bb2HDM_MZp1500_MA0300_match.root";
    //string outputRootFile=Form("nTuple/%s_boost_match.root",a0mass.Data());
    
    TFile* outFile = new TFile(outputRootFile.data(),"recreate");
    h_hM->Write();
    h_a0M->Write();
    h_hpt->Write();
    h_a0pt->Write();
    h_htau->Write();
    h_a0tau->Write();
    h_heta->Write();
    h_a0eta->Write();
    outFile->Close();
    
}
