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
#define doGenMatch      true
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
bool sortListbyPt(vector<float> a, vector<float> b) {return a[3]>b[3];}
bool sortListbyM(vector<float> a, vector<float> b) {return abs(a[2])<abs(b[2]);}
bool sortListbyPtZp(vector<float> a, vector<float> b) {return a[4]>b[4];}
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
void anbb2HDM_genmatch_4jet(int w=0, std::string inputFile="2HDM_MZp1000_MA0300_re.root"){
    
    setNCUStyle(true);
    //gStyle->SetOptTitle(0);
    //gStyle->SetOptStat(0001111101.);
    gStyle->SetOptStat(0);

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
    isBG = false;
    TCanvas* c1 = new TCanvas("c1","",600,600);
    TH1F* h_allEvent = new TH1F("h_allEvent","h_allEvent",10,-0.5,9);
    
    
    TH1F* h_hNcandi = new TH1F("h_hNcandi", "h_higgs_NcombinationJets", 10,-0.5,9.5);
    TH1F* h_a0Ncandi = new TH1F("h_a0Ncandi", "h_A0_NcombinationJets", 10,-0.5,9.5);
    TH1F* h_zpNcandi = new TH1F("h_zpNcandi", "h_Zp_NcombinationJets", 10,-0.5,9.5);

    //TH1F* h_hM = new TH1F("h_higgsM", "h_higgsM_reco", 75,100,175);
    TH1F* h_hM = new TH1F("h_higgsM", "h_higgsM_reco", 75,-100,50);
    TH1F* h_hM_ori = new TH1F("h_higgsM_ori", "h_higgsM_reco", 75,100,175);
    TH1F* h_hM_all = new TH1F("h_higgsM_all", "h_higgsM_all", 100,50,150);
    TH1F* h_a0M = new TH1F("h_a0M", "h_A0M_reco", 40,-300,+100);
    //TH1F* h_a0M = new TH1F("h_a0M", "h_A0M_reco", 24,A0mass.Atof()-60,A0mass.Atof()+60);
    TH1F* h_a0M_all = new TH1F("h_A0M_all", "h_A0M_all", 100,0,500);
    TH1F* h_zpM = new TH1F("h_ZpM", "h_ZpM_reco", 60,Zpmass.Atof()-400,Zpmass.Atof()+200);
    TH1F* h_zpM_ori = new TH1F("h_ZpM_ori", "h_ZpM_reco", 24,Zpmass.Atof()-120,Zpmass.Atof()+120);
    
    Int_t nPass[20]={0};
    int maxHIndexNum = 0, maxZpIndexNum = 0;
    
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        
        if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        data.GetEntry(jEntry);
        //if (jEntry>50) break; 
        nPass[0]++;
        
        //0. has a good vertex
        int nGenJet =  data.GetInt("THINnJet");
        if(nGenJet<4) continue;
        nPass[1]++;
        
        TClonesArray* genjetP4 =  (TClonesArray*) data.GetPtrTObject("THINjetP4");
        
        const int nCandidates = 100;
        vector<vector<float>> HindexList, A0indexList, ZpindexList;
        vector<TLorentzVector*> genHA0Par;
        vector<float> row(5,-99);
        ZpindexList.assign(nCandidates,row);
        
        if (!isBG||doGenMatch) {
            TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
            for (int i=0;i<4;i++) genHA0Par.push_back((TLorentzVector*)genParP4->At(genPar[jEntry][i]));
        }
        
        // reco higgs
        float genMatchHiggs[3]={0}, genMatchA0[3] = {0};
        float diJetM = 0, mindiJetM = 0;
        for(int ij=0; ij < nGenJet; ij++) {
            TLorentzVector* thisJet = (TLorentzVector*)genjetP4->At(ij);
            for (int jj=0;jj<ij;jj++) {
                TLorentzVector* thatJet = (TLorentzVector*)genjetP4->At(jj);
                if (!isBG||doGenMatch) {
                    bool genParA=false, genParB=false;
                    if (genHA0Par[0]->DeltaR(*thisJet)<0.4&&genHA0Par[1]->DeltaR(*thatJet)<0.4) genParA=true;
                    if (genHA0Par[1]->DeltaR(*thisJet)<0.4&&genHA0Par[0]->DeltaR(*thatJet)<0.4) genParB=true; 
                    if (!(genParA||genParB)) continue;
                }

                diJetM = (*thisJet+*thatJet).M();
                if (abs(genMatchHiggs[2])>abs(diJetM-125)) {
                    genMatchHiggs[0] = ij;
                    genMatchHiggs[1] = jj;
                    genMatchHiggs[2] = diJetM-125;

                }
                //h_hM_all->Fill(diJetM);
                HindexList.push_back({(float)ij,(float)jj,(float)(*thisJet+*thatJet).M()-125,(float)(*thisJet+*thatJet).Pt()});
            }
        } // end of outer loop jet
        nPass[2]++;
        sort(HindexList.begin(),HindexList.end(),sortListbyPt);
        h_hNcandi->Fill(HindexList.size());
        
        
        // reco A0
        mindiJetM = 0;
        diJetM = 0;
        for(int ij=0; ij < nGenJet; ij++) {
            
            TLorentzVector* thisJet = (TLorentzVector*)genjetP4->At(ij);
            for (int jj=0;jj<ij;jj++) {
                TLorentzVector* thatJet = (TLorentzVector*)genjetP4->At(jj);
                if (!isBG||doGenMatch) {
                    bool genParA=false, genParB=false;
                    if (genHA0Par[2]->DeltaR(*thisJet)<0.4&&genHA0Par[3]->DeltaR(*thatJet)<0.4) genParA=true;
                    if (genHA0Par[3]->DeltaR(*thisJet)<0.4&&genHA0Par[2]->DeltaR(*thatJet)<0.4) genParB=true; 
                    if (!isBG && !(genParA||genParB)) continue;
                }
                diJetM = (*thisJet+*thatJet).M();
                if (abs(genMatchA0[2])>abs(diJetM-A0mass.Atof())) {
                    genMatchA0[0] = ij;
                    genMatchA0[1] = jj;
                    genMatchA0[2] = diJetM-A0mass.Atof();

                }
                A0indexList.push_back({(float)ij,(float)jj,(float)((*thisJet+*thatJet).M()-A0mass.Atof()),(float)((*thisJet+*thatJet).Pt())});
            }
        } // end of outer loop jet
        sort(A0indexList.begin(),A0indexList.end(),sortListbyPt);
        h_a0Ncandi->Fill(A0indexList.size());
        nPass[3]++;
        if (A0indexList.size()!=0&&HindexList.size()==0) {
            h_a0M->Fill(A0indexList[0][2]);
            continue;
        }
        else if (HindexList.size()!=0&&A0indexList.size()==0) {
            h_hM->Fill(HindexList[0][2]);
            continue;
        }
        float zpM = -999;
        bool isSave = false;
        for (int i=0 ;i<HindexList.size();i++) {
            for (int j=0;j<A0indexList.size();j++) {
                bool idCrash = false;
                for (int m=0;m<2;m++) for (int n=0;n<2;n++) if (HindexList[i][m]==A0indexList[j][n]) idCrash = true;
                if (idCrash) continue;
                TLorentzVector* bJet0 = (TLorentzVector*)genjetP4->At(HindexList[i][0]);
                TLorentzVector* bJet1 = (TLorentzVector*)genjetP4->At(HindexList[i][1]);
                TLorentzVector* bJet2 = (TLorentzVector*)genjetP4->At(A0indexList[j][0]);
                TLorentzVector* bJet3 = (TLorentzVector*)genjetP4->At(A0indexList[j][1]);
                zpM = (*bJet0+*bJet1+*bJet2+*bJet3).M();
                h_hM->Fill(HindexList[i][2]);
                h_a0M->Fill(A0indexList[j][2]);
                h_zpM->Fill(zpM);
                isSave = true;
                if (isSave) break;
            }
            if (isSave) break;   
        }
    } // end of loop over entries
    
    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    
    for(int i=0;i<20;i++) if(nPass[i]>0) std::cout << "nPass[" << i << "] = " << nPass[i] << std::endl;
    efferr(nPass[4],nTotal);

    if (!isBG) 
    {
        string pdfName;
        if (doGenMatch) pdfName = Form("anGenMatchJet_bb2HDM_MZp%s_MA0%s_width.pdf",Zpmass.Data(),A0mass.Data());
        c1->Print((pdfName+"[").data());
        h_zpM->Draw("hist");
        c1->Print(pdfName.data());
        //h_hM_ori->Draw("hist");
        h_hM->Draw("hist");
        c1->Print(pdfName.data());
        //h_hM_all->Draw("hist");
        //c1->Print(pdfName.data());
        h_a0M->Draw("hist");
        c1->Print(pdfName.data());
        //h_a0M_all->Draw("hist");
        //c1->Print(pdfName.data());
        c1->Print((pdfName+"]").data());
    }
    
    string fileName = "fitwidth.root"; 
    if (!saveEffi||true) {
        TFile* outputFile = new TFile(fileName.data(),"recreate");
        h_zpM->Write();
        h_hM->Write();
        h_a0M->Write();
        
        outputFile->Close();
    }
    /*
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
