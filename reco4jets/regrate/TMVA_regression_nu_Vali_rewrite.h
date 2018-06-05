//#include "Utilities.h"
#include "TMVA/Reader.h"


struct TMVAinputInfo{

	Int_t    nJet;
	Int_t    nVtx;
	
	Float_t  pfMET;
	Float_t  pfMETPhi;

	Int_t*   jetNCH;
	Float_t* jetMt;
	Float_t* jetNHF;
	Float_t* jetCHF;
	Float_t* jetNEF;
	Float_t* jetVtxPt;
	Float_t* jetVtxMass;
	Float_t* jetVtx3DVal;
	Float_t* jetVtxNtrks;
	Float_t* jetVtx3DSig;
	Float_t* CISVV2;

	TLorentzVector *rawjet[35];
	TLorentzVector *ak4jet[35];
	TLorentzVector *leadjet[35];
	TLorentzVector *lepjet[35];
	TLorentzVector *Parjet[30];

public:
	TMVAinputInfo() {}
	TMVAinputInfo(TreeReader &data, bool usePar = false) {
		nJet         = data.GetInt("THINnJet");
		nVtx         = data.GetInt("nVtx");
		jetMt        = data.GetPtrFloat("THINjetMt");
		jetNHF       = data.GetPtrFloat("THINjetNHadEF");
		jetCHF       = data.GetPtrFloat("THINjetCHadEF");
		jetNEF       = data.GetPtrFloat("THINjetNEmEF");
		jetNCH       = data.GetPtrInt("THINjetCMulti");
		jetVtxPt     = data.GetPtrFloat("THINjetVtxPt");
		jetVtxMass   = data.GetPtrFloat("THINjetVtxMass");
		jetVtx3DVal  = data.GetPtrFloat("THINjetVtx3DVal");
		jetVtxNtrks  = data.GetPtrFloat("THINjetVtxNtrks");
		jetVtx3DSig  = data.GetPtrFloat("THINjetVtx3DSig");
		pfMET        = data.GetFloat("pfMetCorrPt");
		pfMETPhi     = data.GetFloat("pfMetCorrPhi");
		CISVV2 		 = data.GetPtrFloat("THINjetCISVV2");
		
		static TClonesArray *ak4jetP4, *rawjetP4, *THINjetLeadTrackP4, *LepjetP4, *genParP4;
		ak4jetP4            =  (TClonesArray*) data.GetPtrTObject("THINunCorrJetP4");
		rawjetP4            =  (TClonesArray*) data.GetPtrTObject("THINjetP4");
		THINjetLeadTrackP4  =  (TClonesArray*) data.GetPtrTObject("THINjetLeadTrackP4");
		LepjetP4            =  (TClonesArray*) data.GetPtrTObject("THINjetLepTrackP4");
		static int i;
		for (i=0;i<nJet;i++) {
			ak4jet[i]       =  (TLorentzVector*)ak4jetP4->At(i);
			rawjet[i]       =  (TLorentzVector*)rawjetP4->At(i);
			leadjet[i]      =  (TLorentzVector*)THINjetLeadTrackP4->At(i);
			lepjet[i]       =  (TLorentzVector*)LepjetP4->At(i);
		}
		if (usePar) {
			genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
			for (i=0;i<30;i++) Parjet[i] = (TLorentzVector*)genParP4->At(i);
		}
	}
};
//2017_3_1 update==============
#ifndef TMVA_REGRESSION
#define TMVA_REGRESSION
inline float deltaR(float eta1,float phi1, float eta2,float phi2) {
	return sqrt( (eta2-eta1)*(eta2-eta1) + (phi1-phi2)*(phi1-phi2) );
}
inline float deltaPhi(float phi1, float phi2) {
	return (abs(phi2-phi1)<(Float_t)TMath::Pi())? abs(phi2-phi1): 2*TMath::Pi()-abs(phi2-phi1);
}
inline float calChiSquare(float Mh, float MA0, float MZp) {
    //return TMath::Power((Mh-125)/7.6,2) + TMath::Power(((MA0-300)/22),2) + TMath::Power((MZp-1000)/41,2);
    return TMath::Power((Mh-125)/14.48,2) + TMath::Power(((MA0-300)/37.86),2) + TMath::Power((MZp-1000)/83.62,2);
}
#endif
// structure, new one
// LL


Float_t TMVA_15plus3_jetGenJet_nu_3_1(TMVAinputInfo &info, Int_t i, Float_t jjDR,int mode = -1) 
{
	if (i>=info.nJet) {
		cout << "get out of range index" << endl;
		return 1.;
	}

  	// classification variables
	static float Jet_pt_, Jet_corr_,  Jet_eta_, Jet_mt_;
	static float Jet_leadTrackPt_, Jet_leptonPtRel_, Jet_leptonPt_, Jet_leptonDeltaR_;
	static float Jet_neHEF_, Jet_neEmEF_, Jet_chMult_, Jet_vtxPt_, Jet_vtxMass_, Jet_vtx3dL_, Jet_vtxNtrk_, Jet_vtx3deL_;
	static float Jet_PFMET_, Jet_METDPhi_, Jet_JetDR_;
	static float nVtx_;
  	// MVA classifiers for b-jet
	static TMVA::Reader* tmvaReader = NULL;
	static TMVA::Reader* tmvaReader_L = NULL;
	static TMVA::Reader* tmvaReader_T = NULL;
  	// one-time MVA initialization
	if (mode == 1 && !tmvaReader_L) {
		
		tmvaReader_L = new TMVA::Reader("!Color:Silent");
		tmvaReader_L->AddVariable("Jet_pt",             	&Jet_pt_);
		tmvaReader_L->AddVariable("Jet_eta",            	&Jet_eta_);
		tmvaReader_L->AddVariable("Jet_mt",             	&Jet_mt_);
		tmvaReader_L->AddVariable("Jet_leadTrackPt",    	&Jet_leadTrackPt_);
		tmvaReader_L->AddVariable("Jet_leptonPtRel_new",  	&Jet_leptonPtRel_);
		tmvaReader_L->AddVariable("Jet_leptonPt",       	&Jet_leptonPt_);
		tmvaReader_L->AddVariable("Jet_leptonDeltaR",   	&Jet_leptonDeltaR_);
		tmvaReader_L->AddVariable("Jet_neHEF",				&Jet_neHEF_);
		tmvaReader_L->AddVariable("Jet_neEmEF",         	&Jet_neEmEF_);
		tmvaReader_L->AddVariable("Jet_vtxPt",          	&Jet_vtxPt_);
		tmvaReader_L->AddVariable("Jet_vtxMass",        	&Jet_vtxMass_);
		tmvaReader_L->AddVariable("Jet_vtx3dL",         	&Jet_vtx3dL_);
		tmvaReader_L->AddVariable("Jet_vtxNtrk",        	&Jet_vtxNtrk_);
		tmvaReader_L->AddVariable("Jet_vtx3deL_new",      	&Jet_vtx3deL_);
		tmvaReader_L->AddVariable("nGoodPVs",            	&nVtx_);
		tmvaReader_L->AddVariable("Jet_PFMET",            	&Jet_PFMET_);
		tmvaReader_L->AddVariable("Jet_METDPhi",          	&Jet_METDPhi_);
		tmvaReader_L->AddVariable("Jet_JetDR",          	&Jet_JetDR_);
    	// read weight files
		tmvaReader_L->BookMVA("BDTG method L","BDTG_15plus3_jetGenJet_nu_leading_summer16_2_26.weights.xml");
	}
	else if (mode == 0 && !tmvaReader_T) {
		tmvaReader_T = new TMVA::Reader("!Color:Silent");
		tmvaReader_T->AddVariable("Jet_pt",             	&Jet_pt_);
		tmvaReader_T->AddVariable("Jet_eta",            	&Jet_eta_);
		tmvaReader_T->AddVariable("Jet_mt",             	&Jet_mt_);
		tmvaReader_T->AddVariable("Jet_leadTrackPt",    	&Jet_leadTrackPt_);
		tmvaReader_T->AddVariable("Jet_leptonPtRel_new",  	&Jet_leptonPtRel_);
		tmvaReader_T->AddVariable("Jet_leptonPt",       	&Jet_leptonPt_);
		tmvaReader_T->AddVariable("Jet_leptonDeltaR",   	&Jet_leptonDeltaR_);
		tmvaReader_T->AddVariable("Jet_neHEF",				&Jet_neHEF_);
		tmvaReader_T->AddVariable("Jet_neEmEF",         	&Jet_neEmEF_);
		tmvaReader_T->AddVariable("Jet_vtxPt",          	&Jet_vtxPt_);
		tmvaReader_T->AddVariable("Jet_vtxMass",        	&Jet_vtxMass_);
		tmvaReader_T->AddVariable("Jet_vtx3dL",         	&Jet_vtx3dL_);
		tmvaReader_T->AddVariable("Jet_vtxNtrk",        	&Jet_vtxNtrk_);
		tmvaReader_T->AddVariable("Jet_vtx3deL_new",      	&Jet_vtx3deL_);
		tmvaReader_T->AddVariable("nGoodPVs",            	&nVtx_);
		tmvaReader_T->AddVariable("Jet_PFMET",            	&Jet_PFMET_);
		tmvaReader_T->AddVariable("Jet_METDPhi",          	&Jet_METDPhi_);
		tmvaReader_T->AddVariable("Jet_JetDR",          	&Jet_JetDR_);
    	// read weight files
		tmvaReader_T->BookMVA("BDTG method T","BDTG_15plus3_jetGenJet_nu_trailing_summer16_2_26.weights.xml");
 	} // one-time initialization
	else if (!tmvaReader) {
		tmvaReader = new TMVA::Reader("!Color:Silent");
		tmvaReader->AddVariable("Jet_pt",             	&Jet_pt_);
		tmvaReader->AddVariable("Jet_eta",            	&Jet_eta_);
		tmvaReader->AddVariable("Jet_mt",             	&Jet_mt_);
		tmvaReader->AddVariable("Jet_leadTrackPt",    	&Jet_leadTrackPt_);
		tmvaReader->AddVariable("Jet_leptonPtRel_new",  &Jet_leptonPtRel_);
		tmvaReader->AddVariable("Jet_leptonPt",       	&Jet_leptonPt_);
		tmvaReader->AddVariable("Jet_leptonDeltaR",   	&Jet_leptonDeltaR_);
		tmvaReader->AddVariable("Jet_neHEF",			&Jet_neHEF_);
		tmvaReader->AddVariable("Jet_neEmEF",         	&Jet_neEmEF_);
		tmvaReader->AddVariable("Jet_vtxPt",          	&Jet_vtxPt_);
		tmvaReader->AddVariable("Jet_vtxMass",        	&Jet_vtxMass_);
		tmvaReader->AddVariable("Jet_vtx3dL",         	&Jet_vtx3dL_);
		tmvaReader->AddVariable("Jet_vtxNtrk",        	&Jet_vtxNtrk_);
		tmvaReader->AddVariable("Jet_vtx3deL_new",      &Jet_vtx3deL_);
		tmvaReader->AddVariable("nGoodPVs",            	&nVtx_);
		tmvaReader->AddVariable("Jet_PFMET",            &Jet_PFMET_);
		tmvaReader->AddVariable("Jet_METDPhi",          &Jet_METDPhi_);
		tmvaReader->AddVariable("Jet_JetDR",          	&Jet_JetDR_);
    	// read weight files
		tmvaReader->BookMVA("BDTG method","BDTG_15plus3_jetGenJet_nu_summer16_3_1.weights.xml");
	}
  	// set MVA variables
  	Jet_pt_           = info.ak4jet[i]->Pt();//Jet pT
  	Jet_corr_         = info.ak4jet[i]->Pt()/info.rawjet[i]->Pt();//JEC
  	nVtx_             = info.nVtx;//nPVs
  	Jet_eta_          = info.ak4jet[i]->Eta();//Jet ¦Ç
  	Jet_mt_           = info.jetMt[i];//Jet transverse mass
  	Jet_leadTrackPt_  = (info.leadjet[i]->Pt() < 0) ? 0. : info.leadjet[i]->Pt();//pTLeadTrk, transverse momentum of the leading track in the jet
  	static TLorentzVector TJet, TJetLep;
  	//TJet.SetPtEtaPhiE(info.ak4jet[i]->Pt(),info.ak4jet[i]->Eta(),info.ak4jet[i]->Phi(),info.ak4jet[i]->E());
  	//TJetLep.SetPtEtaPhiE(jetLepTrackPt[i],jetLepTrackEta[i],jetLepTrackPhi[i],jetLepTrackPt[i]*cosh(jetLepTrackEta[i]));
  	TVector3 TJetAxis(info.ak4jet[i]->Vect().X(),info.ak4jet[i]->Vect().Y(),info.ak4jet[i]->Vect().Z());
  	Jet_leptonPtRel_  = (info.lepjet[i]->Pt() > 0 && info.lepjet[i]->Pt() < 99900) ?  info.lepjet[i]->Perp(TJetAxis)	:0.;//Soft Lepton pTRel, relative transverse momentum of soft lepton candidate in the jet;
  	Jet_leptonPt_     = (info.lepjet[i]->Pt() < 0 || info.lepjet[i]->Pt() > 99900) ? 0. : info.lepjet[i]->Pt();//Soft Lepton pT, transverse momentum of soft lepton candidate in the jet;
  	Jet_leptonDeltaR_ = (info.lepjet[i]->Pt() < 0 || info.lepjet[i]->Pt() > 99900) ? 0. : info.lepjet[i]->DeltaR(*info.ak4jet[i]);//Soft Lepton dR, relative distance in the ¦Ç-phi space of soft lepton candidate in the jet and the jet;
  	Jet_neHEF_        = (info.jetNHF[i] < 1.) ? info.jetNHF[i]:1.;//Neutral hadron energy fraction
  	Jet_neEmEF_       = (info.jetNEF[i] < 1.) ? info.jetNEF[i]:1.;//Photon energy fraction
  	//Jet_chMult_       = (Float_t) jetNCH[i];//total number of jet constituents
  	Jet_vtxPt_        = info.jetVtxPt[i];//SecVtxPt, pT of the jet secondary vertex
  	Jet_vtxMass_      = info.jetVtxMass[i];//SecVtxM, Mass of the jet secondary vertex
  	Jet_vtx3dL_       = (info.jetVtx3DVal[i] > 0) ? info.jetVtx3DVal[i] : 0.;//SecVtxdL, the 3-d flight length of the jet secondary vertex
  	Jet_vtxNtrk_      = (Int_t) info.jetVtxNtrks[i];//SecVtxNtrk, track multiplicity of the reconstructed secondary vertex
  	Jet_vtx3deL_	  = (info.jetVtx3DSig[i] > 0) ? Jet_vtx3dL_/info.jetVtx3DSig[i] : 0.;//SecVtxdeL, Error on the 3-d flight length of the jet secondary vertex
  	Jet_PFMET_        = info.pfMET;
  	Jet_METDPhi_      = deltaPhi(info.pfMETPhi,info.ak4jet[i]->Phi());
  	Jet_JetDR_        = jjDR;
	if (mode == 1) 		return tmvaReader_L->EvaluateRegression("BDTG method L")[0];
	else if (mode == 0) return tmvaReader_T->EvaluateRegression("BDTG method T")[0];
	else 				return tmvaReader->EvaluateRegression("BDTG method")[0];
}


struct HA0JetInfo_base
{
	TLorentzVector* Hjet1   = NULL;
	TLorentzVector* Hjet2   = NULL;
	TLorentzVector* A0jet1  = NULL;
	TLorentzVector* A0jet2  = NULL;
	bool dodelete = false;
public:
	void remove(TLorentzVector* v1) {
		if (v1) {
			v1 = NULL;
		}
	}
	HA0JetInfo_base(bool b=false) {
		dodelete = b;
	}
	~HA0JetInfo_base(){
		if (dodelete&&false) {
			remove(Hjet1);
			remove(Hjet2);
			remove(A0jet1);
			remove(A0jet2);
		}
	}
};

struct HA0JetInfo {
	HA0JetInfo_base Jets;
	HA0JetInfo_base weightJets{true};
	HA0JetInfo_base weightJetsLT{true};

private:
	bool isCopy = false;
	float weight[4] 	= {-99,-99,-99,-99};
	float weightLT[4] 	= {-99,-99,-99,-99};
	int   ind[4]    	= {-1,-1,-1,-1};
	float vCISVV2[4]	= {-1,-1,-1,-1};
private:
	Float_t CISVV2CUT_L = 0.5426;
    Float_t CISVV2CUT_M = 0.8484;
    Float_t CISVV2CUT_T = 0.9535;
public:

	TLorentzVector* applyWeight(TLorentzVector *refv1, float weight) {
		if (!refv1) {
			cout << "give null" << endl;
			return NULL;
		}
		TLorentzVector* v1 = new TLorentzVector();
		v1->SetPtEtaPhiE(refv1->Pt()*weight,refv1->Eta(),refv1->Phi(),refv1->E()*weight);
		//v1->SetPtEtaPhiM(refv1->Pt()*weight,refv1->Eta(),refv1->Phi(),refv1->M());
		return v1;
	}
	HA0JetInfo() {
	}
	HA0JetInfo(TMVAinputInfo &info, int hi1, int hi2, int a0i1, int a0i2) {
		
		if (info.ak4jet[hi1]->Pt() < info.ak4jet[hi2]->Pt()) swap(hi1,hi2);
		if (info.ak4jet[a0i1]->Pt() < info.ak4jet[a0i2]->Pt()) swap(a0i1,a0i2);
		ind[0] = hi1;	// leading
		ind[1] = hi2;
		ind[2] = a0i1;	// leading
		ind[3] = a0i2;
		Jets.Hjet1 = info.ak4jet[hi1];
		Jets.Hjet2 = info.ak4jet[hi2];
		Jets.A0jet1 = info.ak4jet[a0i1];
		Jets.A0jet2 = info.ak4jet[a0i2];
		// gen Par
		/*
		Jets.Hjet1 = info.Parjet[hi1];
		Jets.Hjet2 = info.Parjet[hi2];
		Jets.A0jet1 = info.Parjet[a0i1];
		Jets.A0jet2 = info.Parjet[a0i2];
		*/
		// calculate weight
		if (weight[0]<0) weight[0] = TMVA_15plus3_jetGenJet_nu_3_1(info, hi1,  Jets.Hjet1->DeltaR(*Jets.Hjet2));
		if (weight[1]<0) weight[1] = TMVA_15plus3_jetGenJet_nu_3_1(info, hi2,  Jets.Hjet1->DeltaR(*Jets.Hjet2));
		if (weight[2]<0) weight[2] = TMVA_15plus3_jetGenJet_nu_3_1(info, a0i1, Jets.A0jet1->DeltaR(*Jets.A0jet2));
		if (weight[3]<0) weight[3] = TMVA_15plus3_jetGenJet_nu_3_1(info, a0i2, Jets.A0jet1->DeltaR(*Jets.A0jet2));
		
		if (weightLT[0]<0) weightLT[0] = TMVA_15plus3_jetGenJet_nu_3_1(info, hi1,  Jets.Hjet1->DeltaR(*Jets.Hjet2),		1);
		if (weightLT[1]<0) weightLT[1] = TMVA_15plus3_jetGenJet_nu_3_1(info, hi2,  Jets.Hjet1->DeltaR(*Jets.Hjet2), 	0);
		if (weightLT[2]<0) weightLT[2] = TMVA_15plus3_jetGenJet_nu_3_1(info, a0i1, Jets.A0jet1->DeltaR(*Jets.A0jet2),	1);
		if (weightLT[3]<0) weightLT[3] = TMVA_15plus3_jetGenJet_nu_3_1(info, a0i2, Jets.A0jet1->DeltaR(*Jets.A0jet2),	0);
		
		// weight jets
		weightJets.dodelete = true;
		weightJets.Hjet1 = 	applyWeight(Jets.Hjet1,	weight[0]);
		weightJets.Hjet2 = 	applyWeight(Jets.Hjet2,	weight[1]);
		weightJets.A0jet1 = applyWeight(Jets.A0jet1,weight[2]);
		weightJets.A0jet2 =	applyWeight(Jets.A0jet2,weight[3]);
		weightJetsLT.dodelete = true;
		weightJetsLT.Hjet1 = 	applyWeight(Jets.Hjet1,		weightLT[0]);
		weightJetsLT.Hjet2 = 	applyWeight(Jets.Hjet2,		weightLT[1]);
		weightJetsLT.A0jet1 = 	applyWeight(Jets.A0jet1,	weightLT[2]);
		weightJetsLT.A0jet2 = 	applyWeight(Jets.A0jet2,	weightLT[3]);
		// swap index by weight pt
		if (weightJetsLT.Hjet1->Pt()<weightJetsLT.Hjet2->Pt()) {
			swap(weightJetsLT.Hjet1,weightJetsLT.Hjet2);
			swap(weightLT[0],weightLT[1]);
			swap(ind[0],ind[1]);
		}
		if (weightJetsLT.A0jet1->Pt()<weightJetsLT.A0jet2->Pt()) {
			swap(weightJetsLT.A0jet1,weightJetsLT.A0jet2);
			swap(weightLT[2],weightLT[3]);
			swap(ind[2],ind[3]);
		}
		// CISVV2
		for (int i=0;i<4;i++) if (ind[i]>=0) vCISVV2[i] = info.CISVV2[ind[i]];
	}

	void remove(TLorentzVector* v1) {
		if (v1) {
			delete v1;
			v1 = NULL;
		}
	}
	~HA0JetInfo(){
	}
	void Release(){
		if (weightJets.Hjet1)	remove(weightJets.Hjet1);
		if (weightJets.Hjet2)	remove(weightJets.Hjet2);
		if (weightJets.A0jet1)	remove(weightJets.A0jet1);
		if (weightJets.A0jet2)	remove(weightJets.A0jet2); 
		if (weightJetsLT.Hjet1)	remove(weightJetsLT.Hjet1);
		if (weightJetsLT.Hjet2)	remove(weightJetsLT.Hjet2);
		if (weightJetsLT.A0jet1)	remove(weightJetsLT.A0jet1);
		if (weightJetsLT.A0jet2)	remove(weightJetsLT.A0jet2);
	}
	inline float Mh(int H = 0) 	{
		TLorentzVector Hdijet;
		if (H==1) Hdijet =  *(weightJets.Hjet1) + *(weightJets.Hjet2);
		else if (H==2) Hdijet = *(weightJetsLT.Hjet1) + *(weightJetsLT.Hjet2);
		else Hdijet = *(Jets.Hjet1) + *(Jets.Hjet2);
		return Hdijet.M();
	}
	inline float MA0(int A0 = 0,bool print=false)    	{
		TLorentzVector A0dijet;
		if (print) cout << "test: " << weightJetsLT.A0jet1 << "\t" << weightJetsLT.A0jet2 << endl;
		if (A0==1) {
			A0dijet =  *(weightJets.A0jet1) + *(weightJets.A0jet2);
		}
		else if (A0==2) {
			A0dijet = *(weightJetsLT.A0jet1) + *(weightJetsLT.A0jet2);
		}
		else {
			A0dijet = *(Jets.A0jet1) + *(Jets.A0jet2);
		}

		return A0dijet.M();
	}
	inline float MZp(int H = 0, int A0 = 0)    	{
		TLorentzVector Hdijet, A0dijet;
		if (H==1) Hdijet =  *weightJets.Hjet1 + *weightJets.Hjet2;
		else if (H==2) Hdijet = *weightJetsLT.Hjet1 + *weightJetsLT.Hjet2;
		else Hdijet = *Jets.Hjet1 + *Jets.Hjet2;
		if (A0==1) A0dijet =  *weightJets.A0jet1 + *weightJets.A0jet2;
		else if (A0==2) A0dijet = *weightJetsLT.A0jet1 + *weightJetsLT.A0jet2;
		else A0dijet = *Jets.A0jet1 + *Jets.A0jet2;
		return (Hdijet+A0dijet).M();
	}
	inline float xSqure(int H = 0, int A0 = 0) 	{
		float x2 = calChiSquare(Mh(H),MA0(A0),MZp(H,A0));
		return x2;
	}
	inline float xSquare_Xpar(int i=0,float pow=1) {
		//if (i==0) cout << Form("%f",Mh(2)) << endl;
		//else if (i==1) cout << Form("%f",TMath::Power(((MA0(2)-300)/37.86),2)) << endl;
		//else if (i==2) cout << Form("%f",TMath::Power((MZp(2)-1000)/83.62,2)) << endl;
		//else cout << "Nan" << endl;
		if (i==0) return TMath::Power((Mh(2)-125)/14.48,pow);
		else if (i==1) return TMath::Power(((MA0(2)-300)/37.86),pow);
		else if (i==2) return TMath::Power((MZp(2)-1000)/83.62,pow);
		else return -1;
	}
	inline float  getWeight(int i) 	{return (i<4&&i>=0)?weight[i]:-1;}
	inline float  getWeightLT(int i) 	{return (i<4&&i>=0)?weightLT[i]:-1;}
	inline int* getIndex()			{return ind;}
	inline int  getIndex(int i)  	{return (i<4&&i>=0)?ind[i]:-1;}

	inline float* getCISVV2()		{return vCISVV2;}
	inline float  getCISVV2(int i)	{return (i<4&&i>=0)?vCISVV2[i]:-1;}
	inline bool isCISVV2_L(int i)	{return getCISVV2(i)>CISVV2CUT_L;}
	inline bool isCISVV2_M(int i)	{return getCISVV2(i)>CISVV2CUT_M;}
	inline bool isCISVV2_T(int i)	{return getCISVV2(i)>CISVV2CUT_T;}

	inline float getMinCISVV2(int i) 	{
		if (i==0) return min(getCISVV2(0),getCISVV2(1));
		else if (i==1) return min(getCISVV2(2),getCISVV2(3));
		else return -1;
	}
	inline bool isMinCISVV2_L(int i) {return getMinCISVV2(i)>CISVV2CUT_L;} 
	inline bool isMinCISVV2_M(int i) {return getMinCISVV2(i)>CISVV2CUT_M;} 
	inline bool isMinCISVV2_T(int i) {return getMinCISVV2(i)>CISVV2CUT_T;} 
	
	inline float getMeanCISVV2(int i) {
		if (i==0) return (getCISVV2(0)+getCISVV2(1))/2;
		else if (i==1) return (getCISVV2(2)+getCISVV2(3))/2;
		else return -1;
	}
	inline bool isMeanCISVV2_L(int i) {return getMeanCISVV2(i)>CISVV2CUT_L;}
	inline bool isMeanCISVV2_M(int i) {return getMeanCISVV2(i)>CISVV2CUT_M;}
	inline bool isMeanCISVV2_T(int i) {return getMeanCISVV2(i)>CISVV2CUT_T;}
	inline float Pth()	  	{return (*weightJetsLT.Hjet2+*weightJetsLT.Hjet1).Pt();}
	inline float HdR()    	{return Jets.Hjet1->DeltaR(*Jets.Hjet2);}
	inline float HdEta()	{return abs(Jets.Hjet1->Eta()-Jets.Hjet2->Eta());}
	inline float HdPhi()	{return abs(Jets.Hjet1->DeltaPhi(*Jets.Hjet2));}
	inline float HptAs()	{
		return pow(weightJetsLT.Hjet2->Pt()*HdR()/Mh(),2);
	}
	inline float HptSD() 	{return weightJetsLT.Hjet2->Pt()/(weightJetsLT.Hjet2->Pt()+weightJetsLT.Hjet1->Pt());}
	/*
	inline float Pta0()		{return (*A0jet2+*A0jet1).Pt();}
	inline float PtZp()		{return (*Hjet2+*Hjet1+*A0jet2+*A0jet1).Pt();}
	inline float A0dR()   	{return A0jet1->DeltaR(*A0jet2);}
	inline float ZpdR()   	{return (*Hjet2+*Hjet1).DeltaR(*A0jet2+*A0jet1);}
	*/

	
	/*
	inline float* getWeight()      	{return weight;}
	
	inline TLorentzVector* getHjet(int i)	{
		if (i==0) return Hjet1;
		else if (i==1) return Hjet2;
		else return NULL;
	}
	inline TLorentzVector* getA0jet(int i)	{
		if (i==0) return A0jet1;
		else if (i==1) return A0jet2;
		else return NULL;
	}
	inline TLorentzVector* getHA0jet(int i) 	{
		TLorentzVector* vec = new TLorentzVector();
		if (i==0) *vec = *Hjet1+*Hjet2;
		else if (i==1) *vec = *A0jet1+*A0jet2;
		else vec = NULL;
		return vec;
	}
	*/
};



