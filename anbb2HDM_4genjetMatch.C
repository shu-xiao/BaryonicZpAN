
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
bool sortListbyPtZp(vector<float> a, vector<float> b) {return a[4]>b[4];}

using namespace std;
void anbb2HDM_4genjetMatch(int w, std::string inputFile){
    setNCUStyle(true);
    //get TTree from file ...
    TreeReader data(inputFile.data());
    
    bool upeff = false;
    TString Zpmass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*_MZp}; test=${test1%%_MA0*.root}; echo \"${test}\"",inputFile.data()));
    TString A0mass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*MA0}; test=${test1%%.root}; echo \"${test}\"",inputFile.data()));

    TCanvas* c1 = new TCanvas("c1","",600,600);

    TH1F* h_hpt = new TH1F("h_higgsPt", "h_higgsPt_reco", 20,0,1400);
    
    TH1F* h_hNcandi = new TH1F("h_hNcandi", "h_higgs_NcandidatePairJets", 10,-0.5,9.5);
    TH1F* h_a0Ncandi = new TH1F("h_a0Ncandi", "h_A0_NcandidatePairJets", 10,-0.5,9.5);
    TH1F* h_zpNcandi = new TH1F("h_zpNcandi", "h_Zp_NcandidatePairJets", 10,-0.5,9.5);

    TH1F* h_hMGen = new TH1F("h_higgsM_gen", "h_higgsM_genPar", 50,100,150);
    TH1F* h_hM = new TH1F("h_higgsM", "h_higgsM_genJet", 50,100,150);
    TH1F* h_a0M = new TH1F("h_a0M", "h_A0M_genJet", 24,A0mass.Atof()-60,A0mass.Atof()+60);
    TH1F* h_zpM = new TH1F("h_ZpM", "h_ZpM_genJet", 24,Zpmass.Atof()-120,Zpmass.Atof()+120);
    
    TH1F* h_hDeltaR = new TH1F("h_HiggstobbdeltaR", "h_HiggstobbDeltaR_reco", 30,0,6);
    TH1F* h_hDeltaPhi = new TH1F("h_HiggstobbdeltaPhi", "h_HiggstobbDeltaPhi_reco", 32,0,3.2);
    TH1F* h_hDeltaEta = new TH1F("h_HiggstobbdeltaEta", "h_HiggstobbDeltaEta_reco", 30,0,6);
    
    Int_t nPass[20]={0};
    
    //for(Long64_t jEntry=0; jEntry<50 ;jEntry++){
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        
        if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        data.GetEntry(jEntry);
        
        //0. has enough jet
        
        TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
        int* genMomParId= data.GetPtrInt("genMomParId");
        int* genParId= data.GetPtrInt("genParId");
        vector<vector<float>> HindexList, A0indexList, ZpindexList;
        vector<float> row(5,-99);
        int Hindex[2]={-1,-1}, A0index[2]={-1,-1};
        // create and initialize a 30*5 2D vector
        const int nCandidates = 30;
        HindexList.assign(nCandidates,row);
        A0indexList.assign(nCandidates,row);
        ZpindexList.assign(nCandidates,row);
        //int Hindex[2] = {-1,-1}, A0index[2] = {-1,-1};
        // match higgs
        for(int ij=0; ij < 30; ij++) {
            TLorentzVector* thisParP4 = (TLorentzVector*)genParP4->At(ij);
            if (abs(genParId[ij])!=5) continue;
            if(upeff && thisParP4->Pt()<30) continue;
            if(upeff && fabs(thisParP4->Eta())>2.4) continue;
            if (genMomParId[ij]==25 && Hindex[0]<0) Hindex[0] = ij;
            else if (genMomParId[ij]==25 && Hindex[1]<0) Hindex[1] = ij;
            if (genMomParId[ij]==28 && A0index[0]<0) A0index[0] = ij;
            else if (genMomParId[ij]==28 && A0index[1]<0) A0index[1] = ij;
        } // end of outer loop jet
        if (Hindex[1]<0 || Hindex[0]<0 || A0index[1]<0 || A0index[0]<0) continue;
        nPass[0]++;
        
        //cout << "Event: " << jEntry << endl;
        
        TLorentzVector* HbPar0 = (TLorentzVector*)genParP4->At(Hindex[0]);
        TLorentzVector* HbPar1 = (TLorentzVector*)genParP4->At(Hindex[1]);
        h_hMGen->Fill((*HbPar0+*HbPar1).M());
        TLorentzVector* A0bPar0 = (TLorentzVector*)genParP4->At(A0index[0]);
        TLorentzVector* A0bPar1 = (TLorentzVector*)genParP4->At(A0index[1]);
        //h_hM->Fill((*HbJet0+*HbJet1).M());
        
        // match Higgs
        int nGenJet = data.GetInt("ak4nGenJet");
        TClonesArray* genjetP4 = (TClonesArray*) data.GetPtrTObject("ak4GenJetP4");
        for(int ij=0, iList=0, jList=0; ij < nGenJet; ij++) {
            TLorentzVector* thisJet = (TLorentzVector*)genjetP4->At(ij);
            if(upeff && thisJet->Pt()<30)continue;
            if(upeff && fabs(thisJet->Eta())>2.4)continue;
            for (int jj=0; jj<ij;jj++) {
                bool a0Mp=false, hMp=false;
                TLorentzVector* thatJet = (TLorentzVector*)genjetP4->At(jj);
                float diJetM = (*thisJet+*thatJet).M();
                if (diJetM<140&&diJetM>110) hMp=true;
                if (diJetM<(A0mass.Atof()+50) && diJetM>(A0mass.Atof()-50)) a0Mp=true;
                if (thisJet->DeltaR(*HbPar0)<0.8 && thatJet->DeltaR(*HbPar1)<0.8 && hMp) {
                    HindexList[iList][0] = ij;
                    HindexList[iList][1] = jj;
                    HindexList[iList][2] = (*thisJet+*thatJet).Pt();
                    iList++;
                }
                if (thisJet->DeltaR(*A0bPar0)<0.8 && thatJet->DeltaR(*A0bPar1)<0.8 && a0Mp) {
                    A0indexList[jList][0] = ij;
                    A0indexList[jList][1] = jj;
                    A0indexList[jList][2] = (*thisJet+*thatJet).Pt();
                    jList++;
                }
                if (iList>=nCandidates || jList>=nCandidates) {
                    cout << "Index out of setting" << endl;
                    return;
                }
                
            }
        } // end of loop jet
        
        sort(HindexList.begin(),HindexList.end(),sortListbyPt);
        sort(A0indexList.begin(),A0indexList.end(),sortListbyPt);
        if (HindexList[0][0]<0) {h_hNcandi->Fill(0); continue;}
        nPass[1]++;
        // remove useless rows in 2D vector
        //cout << jEntry << endl;
        for (int i=0;i<HindexList.size();i++) {
            if (HindexList[i][0]<0) {
                h_hNcandi->Fill(i);
                HindexList.erase(HindexList.begin()+i,HindexList.end());
                break;
            }
            //cout << "Hjet" << i << " [ " << HindexList[i][0] << " , " << HindexList[i][1] << " , " << HindexList[i][2] << " ]" << endl;
        }
        
        if (A0indexList[0][0]<0) {h_a0Ncandi->Fill(0); continue;}
        nPass[2]++;
        for (int i=0;i<A0indexList.size();i++) {
            if (A0indexList[i][0]<0) {
                h_a0Ncandi->Fill(i);
                A0indexList.erase(A0indexList.begin()+i,A0indexList.end());
                break;
            }
            cout << "A0jet" << i << " [ " << A0indexList[i][0] << " , " << A0indexList[i][1] << " , " << A0indexList[i][2] << " ]" << endl;
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
                if (zpM>(Zpmass.Atof()+100) || zpM<(Zpmass.Atof()-100)) continue;
                ZpindexList[iList][0] = HindexList[i][0];
                ZpindexList[iList][1] = HindexList[i][1];
                ZpindexList[iList][2] = A0indexList[j][0];
                ZpindexList[iList][3] = A0indexList[j][1];
                ZpindexList[iList][4] = (*bJet0+*bJet1+*bJet2+*bJet3).Pt();
                iList++;
                if (iList>=nCandidates) {
                    cout << "Zp index out of setting" << endl;
                    return;
                }
            }
        }
        if (ZpindexList[0][0]<0) {h_zpNcandi->Fill(0); continue;}
        nPass[3]++;
        sort(ZpindexList.begin(),ZpindexList.end(),sortListbyPtZp);
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
    string pdfName = Form("genMatch_bb2HDM_MZp%s_MA0%s.pdf",Zpmass.Data(),A0mass.Data());
    string pdfNameI = Form("genMatch_bb2HDM_MZp%s_MA0%s.pdf(",Zpmass.Data(),A0mass.Data());
    string pdfNameF = Form("genMatch_bb2HDM_MZp%s_MA0%s.pdf)",Zpmass.Data(),A0mass.Data());
    h_zpM->Draw("hist");
    c1->Print(pdfNameI.data());
    h_hM->Draw("hist");
    c1->Print(pdfName.data());
    h_a0M->Draw("hist");
    c1->Print(pdfName.data());
    h_hMGen->Draw("hist");
    c1->Print(pdfName.data());
    h_hNcandi->Draw("hist");
    c1->Print(pdfName.data());
    h_a0Ncandi->Draw("hist");
    c1->Print(pdfName.data());
    h_zpNcandi->Draw("hist");
    c1->Print(pdfNameF.data());
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
