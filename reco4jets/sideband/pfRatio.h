
float calChiSquare(float mh, float ma0, float mzp);
struct fourJetInfo;

void efferr(float nsig,float ntotal,float factor=1)
{
    float eff = nsig/ntotal;
    float err = sqrt( (1-eff)*eff/ntotal);
    cout << "efficiency = " << eff*factor << " +- " << err*factor << endl;
}
float caldePhi(float phi1, float phi2) {
    float dePhi = 0;
    if (abs(phi1-phi2)>TMath::Pi()) dePhi = 2*TMath::Pi() - abs(phi1-phi2);
    else dePhi = abs(phi1-phi2);
    return dePhi;
}
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
bool sortListbyPt(vector<float> a, vector<float> b) {return a[2]>b[2];}
float calChiSquare(float Mh, float MA0, float MZp) {
    return TMath::Power((Mh-125)/7.6,2) + TMath::Power(((MA0-300)/22),2) + TMath::Power((MZp-1000)/41,2);
}

struct fourJetInfo{
    float Mh;
    float MA0;
    float MZp;
    float Pth;
    float PtA0;
    float PtZp;
    float ChiSquare;
    float CISVV2_HA0[4];
    float minCISVV2[2];
    fourJetInfo() {};
    fourJetInfo(TLorentzVector* hj1,TLorentzVector* hj2, TLorentzVector* a0j1, TLorentzVector* a0j2) {
        Mh   = (*hj1+*hj2).M();
        MA0  = (*a0j1+*a0j2).M();
        MZp  = (*hj1+*hj2+*a0j1+*a0j2).M();
        Pth  = (*hj1+*hj2).Pt();
        PtA0 = (*a0j1+*a0j2).Pt();
        PtZp = (*hj1+*hj2+*a0j1+*a0j2).Pt();
        ChiSquare = calChiSquare(Mh,MA0,MZp);
    }
    fourJetInfo(TreeReader &data, int hi1, int hi2, int a0i1, int a0i2) {
        float *CISVV2 = data.GetPtrFloat("THINjetCISVV2");
        TClonesArray* genjetP4 =  (TClonesArray*) data.GetPtrTObject("THINjetP4");
        TLorentzVector* hj1 = (TLorentzVector*)genjetP4->At(hi1);
        TLorentzVector* hj2 = (TLorentzVector*)genjetP4->At(hi2);
        TLorentzVector* a0j1 = (TLorentzVector*)genjetP4->At(a0i1);
        TLorentzVector* a0j2 = (TLorentzVector*)genjetP4->At(a0i2);
        Mh   = (*hj1+*hj2).M();
        MA0  = (*a0j1+*a0j2).M();
        MZp  = (*hj1+*hj2+*a0j1+*a0j2).M();
        Pth  = (*hj1+*hj2).Pt();
        PtA0 = (*a0j1+*a0j2).Pt();
        PtZp = (*hj1+*hj2+*a0j1+*a0j2).Pt();
        ChiSquare = calChiSquare(Mh,MA0,MZp);
        CISVV2_HA0[0] = CISVV2[hi1];
        CISVV2_HA0[1] = CISVV2[hi2];
        CISVV2_HA0[2] = CISVV2[a0i1];
        CISVV2_HA0[3] = CISVV2[a0i2];
        minCISVV2[0] = min(CISVV2_HA0[0],CISVV2_HA0[1]);
        minCISVV2[1] = min(CISVV2_HA0[2],CISVV2_HA0[3]);
    }
};
