
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
using namespace std;
void resolved_matchBZ4bjet(int w, std::string inputFile){
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

    TH1F* h_oriCIS = new TH1F("h_oriCIS", "h_originCISVV2", 40,0,1);
    TH1F* h_hCIS0 = new TH1F("h_higgsCIS0", "h_higgstobCISVV2_0", 40,0,1);
    TH1F* h_hCIS1 = new TH1F("h_higgsCIS1", "h_higgstobCISVV2_1", 40,0,1);
    TH1F* h_zpCIS0 = new TH1F("h_zpCIS0", "h_zptobCISVV2_0", 40,0,1);
    TH1F* h_zpCIS1 = new TH1F("h_zpCIS1", "h_zptobCISVV2_1", 40,0,1);
    TH1F* h_sumCIS = new TH1F("","",40,0,1);
    TH1F* h_npassCIST = new TH1F("h_npassCIST", "h_zphiggspassCIS_tight", 6,-0.5,5.5);
    TH1F* h_npassCISM = new TH1F("h_npassCISM", "h_zphiggspassCIS_middle", 6,-0.5,5.5);
    TH1F* h_npassCISL = new TH1F("h_npassCISL", "h_zphiggspassCIS_loose", 6,-0.5,5.5);
    TH1F* h_hnpassCIST = new TH1F("h_hnpassCIST", "h_higgspassCIS_tight", 4,-0.5,3.5);
    TH1F* h_hnpassCISM = new TH1F("h_hnpassCISM", "h_higgspassCIS_middle", 4,-0.5,3.5);
    TH1F* h_hnpassCISL = new TH1F("h_hnpassCISL", "h_higgspassCIS_loose", 4,-0.5,3.5);
    TH1F* h_zpnpassCIST = new TH1F("h_zpnpassCIST", "h_zppassCIS_tight", 4,-0.5,3.5);
    TH1F* h_zpnpassCISM = new TH1F("h_zpnpassCISM", "h_zppassCIS_middle", 4,-0.5,3.5);
    TH1F* h_zpnpassCISL = new TH1F("h_zpnpassCISL", "h_zppassCIS_loose", 4,-0.5,3.5);
    
    TH1F* h_hM = new TH1F("h_higgsM", "h_higgsM_match", 30,0,300);
    TH1F* h_zpM = new TH1F("h_ZpM", "h_ZpM_match", 25,0,2000);
    TH1F* h_4bM = new TH1F("h_4bM","h_4bM_match", 25,500,2500);
    
    TH1F* h_hDeltaEta = new TH1F("h_HiggstobbdeltaEta", "h_HiggstobbDeltaEta_match", 25,0,6);
    TH1F* h_zpDeltaEta = new TH1F("h_ZptobbDeltaEta", "h_ZptobbDeltaEta_match", 25,0,6);
    TH1F* h_hDeltaPhi = new TH1F("h_HiggstobbdeltaR", "h_HiggstobbDeltaPhi_match", 25,0,3.2);
    TH1F* h_zpDeltaPhi = new TH1F("h_ZptobbDeltaR", "h_ZptobbDeltaPhi_match", 25,0,3.2);
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
            if (thinJetCSV[k]>=0) h_oriCIS->Fill(thinJetCSV[k]);
            TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(k);
            if (abs(thisJet->Eta())>2.4) continue;
            if (thisJet->Pt()<30) continue;
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
        h_4bM->Fill((*matchHtobJet0+*matchHtobJet1+*matchZptobJet0+*matchZptobJet1).M());
        h_htobbpt0->Fill(matchHtobJet0->Pt());
        h_htobbpt1->Fill(matchHtobJet1->Pt());
        h_zptobbpt0->Fill(matchZptobJet0->Pt());
        h_zptobbpt1->Fill(matchZptobJet1->Pt());

        h_hCIS0->Fill(thinJetCSV[matchHIndex[0]]);
        h_hCIS1->Fill(thinJetCSV[matchHIndex[1]]);
        h_zpCIS0->Fill(thinJetCSV[matchZpIndex[0]]);
        h_zpCIS1->Fill(thinJetCSV[matchZpIndex[1]]);

        h_hDeltaEta->Fill(abs(matchHtobJet0->Eta()-matchHtobJet1->Eta()));
        h_hDeltaPhi->Fill(caldePhi(matchHtobJet0->Phi(),matchHtobJet1->Phi()));
        h_hDeltaR->Fill(matchHtobJet0->DeltaR(*matchHtobJet1));
        h_zpDeltaEta->Fill(abs(matchZptobJet0->Eta()-matchZptobJet1->Eta()));
        h_zpDeltaPhi->Fill(caldePhi(matchZptobJet0->Phi(),matchZptobJet1->Phi()));
        h_zpDeltaR->Fill(matchZptobJet0->DeltaR(*matchZptobJet1));

        // N CISVV2 pass
        int hnpassT = 0, hnpassM = 0, hnpassL = 0;
        int zpnpassT = 0, zpnpassM = 0, zpnpassL = 0;
        const float thinCISVV2cutL = 0.5426; // tight cut
        const float thinCISVV2cutM = 0.8484; // mid cut
        const float thinCISVV2cutT = 0.9535; // tight cut
        for (int i=0;i<2;i++) {
            if (thinJetCSV[matchHIndex[i]]>thinCISVV2cutL) hnpassL++;    
            if (thinJetCSV[matchZpIndex[i]]>thinCISVV2cutL) zpnpassL++;    
            if (thinJetCSV[matchHIndex[i]]>thinCISVV2cutM) hnpassM++;    
            if (thinJetCSV[matchZpIndex[i]]>thinCISVV2cutM) zpnpassM++;    
            if (thinJetCSV[matchHIndex[i]]>thinCISVV2cutT) hnpassT++;    
            if (thinJetCSV[matchZpIndex[i]]>thinCISVV2cutT) zpnpassT++;    
        }
        h_npassCISL->Fill(hnpassL+zpnpassL);
        h_npassCISM->Fill(hnpassM+zpnpassM);
        h_npassCIST->Fill(hnpassT+zpnpassT);
        h_hnpassCISL->Fill(hnpassL);
        h_hnpassCISM->Fill(hnpassM);
        h_hnpassCIST->Fill(hnpassT);
        h_zpnpassCISL->Fill(zpnpassL);
        h_zpnpassCISM->Fill(zpnpassM);
        h_zpnpassCIST->Fill(zpnpassT);
    
    } // end of loop over entries
    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    for(int i=0;i<20;i++) if(nPass[i]>0) std::cout << "nPass[" << i << "] = " << nPass[i] << std::endl;
    efferr(nPass[4],nTotal);
   
    h_sumCIS->Add(h_hCIS0);
    h_sumCIS->Add(h_hCIS1);
    h_sumCIS->Add(h_zpCIS0);
    h_sumCIS->Add(h_zpCIS1);
    h_sumCIS->Draw("hist");
    h_sumCIS->SetLineColor(kBlue);

    h_zpM->Draw("hist");
    c1->Print("thinJetgenMatch.pdf(");
    h_hM->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_4bM->Draw("hist");
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
    h_oriCIS->Draw("hist");
    h_sumCIS->Draw("histsame");
    c1->Print("thinJetgenMatch.pdf");
    h_hCIS0->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_hCIS1->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_zpCIS0->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_zpCIS1->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_npassCIST->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_npassCISM->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_npassCISL->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_hnpassCIST->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_hnpassCISM->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_hnpassCISL->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_zpnpassCIST->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_zpnpassCISM->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_zpnpassCISL->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_hDeltaEta->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_hDeltaPhi->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_hDeltaR->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_zpDeltaEta->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_zpDeltaPhi->Draw("hist");
    c1->Print("thinJetgenMatch.pdf");
    h_zpDeltaR->Draw("hist");
    c1->Print("thinJetgenMatch.pdf)");
}
