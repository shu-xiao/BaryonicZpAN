
#include "iostream"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

using namespace std;
const TString bTag[3] = {"T","M","L"};
struct fourJetCom {
    TString sCom;
    double effi;
    double nB = 0;
    double punzi;
    
    void setName(int i, int j, int k) {
        if ((i+j+k)!=4) return;
        sCom = "";
        for (int n=0;n<i;n++) sCom += bTag[0];
        for (int n=0;n<j;n++) sCom += bTag[1];
        for (int n=0;n<k;n++) sCom += bTag[2];
    }
    void setPunzi() {
        punzi = effi/(1+TMath::Sqrt(nB));
    }
    void Print() {cout << sCom << "\teffi: " << effi << "\tpunzi: " << punzi << "\tnB: " << nB <<endl;}
};
bool sortPunzi(fourJetCom s1, fourJetCom s2) {return s1.punzi>s2.punzi;}
float getL(int nEvent, double xs) {
    return nEvent/xs;
}

void TTLL() {
    string HT_list[9] = {"QCD_TMVA_HT50to100","QCD_TMVA_HT100to200","QCD_TMVA_HT200to300","QCD_TMVA_HT300to500","QCD_TMVA_HT500to700","QCD_TMVA_HT700to1000","QCD_TMVA_HT1000to1500","QCD_TMVA_HT1500to2000","QCD_TMVA_HT2000toInf"};
    TFile* file[9], *fSig;
    file[0] = TFile::Open("QCD_HT50to100.root");
    file[1] = TFile::Open("QCD_HT100to200.root");
    file[2] = TFile::Open("QCD_HT200to300.root");
    file[3] = TFile::Open("QCD_HT300to500.root");
    file[4] = TFile::Open("QCD_HT500to700.root");
    file[5] = TFile::Open("QCD_HT700to1000.root");
    file[6] = TFile::Open("QCD_HT1000to1500.root");
    file[7] = TFile::Open("QCD_HT1500to2000.root");
    file[8] = TFile::Open("QCD_HT2000toInf.root");
    fSig = TFile::Open("treeLead4jet_signal.root");
    const double xsHTbeam[9] = {246400000,27990000,1712000,347700,32100,6831,1207,119.9,25.24};
    const double L2016 = 35.9*1000; // 35.9 fb^-1

    fourJetCom btagCom[15];
    int n = 0;
    for (int i=0;i<=4;i++) for(int j=0;i+j<=4;j++) {btagCom[n].setName(i,j,4-i-j);n++;}
    
    TTree* t,*tSig;
    Int_t nPassL, nPassM, nPassT;
    Int_t nSigPassL, nSigPassM, nSigPassT;
    Int_t nLjet, nMjet, nTjet;
    Int_t nEvent;
    Int_t nEList[9]={0};
    long nTotalEvent = 0;
    TH1F *h_all;
    int nPass[15] = {0};

    // calculate sig effi
    tSig = (TTree*)fSig->Get("tree");
    tSig->SetBranchAddress("nPassL",&nSigPassL);
    tSig->SetBranchAddress("nPassM",&nSigPassM);
    tSig->SetBranchAddress("nPassT",&nSigPassT);
    for (int iE=0;iE<tSig->GetEntries();iE++) {
        tSig->GetEntry(iE);
        n = 0; 
        for (int i=0;i<=4;i++) {
            // loop M, L = 4 - T - M
            for(int j=0, k=0;i+j<=4;j++) {
                k = 4-i-j;
                if (nSigPassT>=i&&nSigPassM>=(i+j)) nPass[n]++;
                n++;
            }
        }
    }
    cout << "end of loop Sig" << endl;
    for (int i=0;i<15;i++) cout << btagCom[i].sCom << "\t" <<  nPass[i] << endl;
    for (int i=0;i<15;i++) btagCom[i].effi = (float)nPass[i]/10000.;
    
    // loop QCD background HT
    for (int xs=0;xs<9;xs++) {
        if (xs==0) continue;
        t = (TTree*)file[xs]->Get("tree");
        t->SetBranchAddress("nPassL",&nPassL);
        t->SetBranchAddress("nPassM",&nPassM);
        t->SetBranchAddress("nPassT",&nPassT);
        h_all = (TH1F*)file[xs]->Get("h_allEvent");
        h_all->Sumw2();
        //nEvent = h_all->Integral();
        nEvent = h_all->GetEntries();
        nEList[xs] = nEvent;
        nTotalEvent += (nEvent*L2016/getL(nEvent,xsHTbeam[xs]));
        for (int i=0;i<15;i++) nPass[i] = 0;
        // loop tree
        for (long iE=0;iE<t->GetEntries();iE++) {
            t->GetEntry(iE);
            int n = 0;
            // loop T
            for (int i=0;i<=4;i++) {
                // loop M, L = 4 - T - M
                for(int j=0, k=0;i+j<=4;j++) {
                    k = 4-i-j;
                    if (nPassT>=i&&nPassM>=(i+j)) nPass[n]++;
                    n++;
                }
            }
            for (n=0;n<15;n++) btagCom[n].nB+= (nPass[n]*L2016/getL(nEvent,xsHTbeam[xs]));

        }
    }
    for (int i=0;i<15;i++) btagCom[i].setPunzi();
    sort(btagCom,btagCom+15,sortPunzi);
    for (int i=0;i<15;i++) btagCom[i].Print();
    cout << "total bk num: " << nTotalEvent << endl;
    for (int i=1;i<9;i++) cout << HT_list[i] << "\t" << nEList[i]<<endl; 
}
