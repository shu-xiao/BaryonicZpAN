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
	

	TLorentzVector *rawjet[35];
	TLorentzVector *ak4jet[35];
	TLorentzVector *leadjet[35];
	TLorentzVector *lepjet[35];
	TMVAinputInfo() {}
	TMVAinputInfo(TreeReader &data) {
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
		
		static TClonesArray *ak4jetP4, *rawjetP4, *THINjetLeadTrackP4, *LepjetP4;
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
    return TMath::Power((Mh-125)/7.6,2) + TMath::Power(((MA0-300)/22),2) + TMath::Power((MZp-1000)/41,2);
}
#endif
// structure, new one
Float_t TMVA_15plus3_jetGenJet_nu_3_1(TMVAinputInfo &info, Int_t i, Float_t jjDR) 
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
  	// one-time MVA initialization
	if (!tmvaReader) {
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
 	} // one-time initialization
  	// set MVA variables
  	Jet_pt_           = info.ak4jet[i]->Pt();//Jet pT
  	Jet_corr_         = info.ak4jet[i]->Pt()/info.rawjet[i]->Pt();//JEC
  	nVtx_             = info.nVtx;//nPVs
  	Jet_eta_          = info.ak4jet[i]->Eta();//Jet η
  	Jet_mt_           = info.jetMt[i];//Jet transverse mass
  	Jet_leadTrackPt_  = (info.leadjet[i]->Pt() < 0) ? 0. : info.leadjet[i]->Pt();//pTLeadTrk, transverse momentum of the leading track in the jet
  	static TLorentzVector TJet, TJetLep;
  	//TJet.SetPtEtaPhiE(info.ak4jet[i]->Pt(),info.ak4jet[i]->Eta(),info.ak4jet[i]->Phi(),info.ak4jet[i]->E());
  	//TJetLep.SetPtEtaPhiE(jetLepTrackPt[i],jetLepTrackEta[i],jetLepTrackPhi[i],jetLepTrackPt[i]*cosh(jetLepTrackEta[i]));
  	TVector3 TJetAxis(info.ak4jet[i]->Vect().X(),info.ak4jet[i]->Vect().Y(),info.ak4jet[i]->Vect().Z());
  	Jet_leptonPtRel_  = (info.lepjet[i]->Pt() > 0 && info.lepjet[i]->Pt() < 99900) ?  info.lepjet[i]->Perp(TJetAxis)	:0.;//Soft Lepton pTRel, relative transverse momentum of soft lepton candidate in the jet;
  	Jet_leptonPt_     = (info.lepjet[i]->Pt() < 0 || info.lepjet[i]->Pt() > 99900) ? 0. : info.lepjet[i]->Pt();//Soft Lepton pT, transverse momentum of soft lepton candidate in the jet;
  	Jet_leptonDeltaR_ = (info.lepjet[i]->Pt() < 0 || info.lepjet[i]->Pt() > 99900) ? 0. : info.lepjet[i]->DeltaR(*info.ak4jet[i]);//Soft Lepton dR, relative distance in the η-phi space of soft lepton candidate in the jet and the jet;
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
  	return tmvaReader->EvaluateRegression("BDTG method")[0];
}

// old one --> fix
Float_t TMVA_15plus3_jetGenJet_nu(TreeReader &data, Int_t i, Float_t jjDR) 
{
	static TClonesArray *ak4jetP4, *rawjetP4, *THINjetLeadTrackP4, *LepjetP4;
	static Int_t nJet;
	static TLorentzVector *RawJet, *thisJet, *leadJet, *lepJet;
	
	ak4jetP4            =  (TClonesArray*) data.GetPtrTObject("THINunCorrJetP4");
	rawjetP4            =  (TClonesArray*) data.GetPtrTObject("THINjetP4");
	THINjetLeadTrackP4  =  (TClonesArray*) data.GetPtrTObject("THINjetLeadTrackP4");
	LepjetP4            =  (TClonesArray*) data.GetPtrTObject("THINjetLepTrackP4");
	nJet                = data.GetInt("THINnJet");
	RawJet  = (TLorentzVector*)rawjetP4->At(i);	

	thisJet = (TLorentzVector*)ak4jetP4->At(i);
	leadJet = (TLorentzVector*)THINjetLeadTrackP4->At(i);
	lepJet = (TLorentzVector*)LepjetP4->At(i);
	
	static Float_t jetPt, jetEn, jetEta, jetPhi, jetRawPt, jetLeadTrackPt, jetLepTrackPt, jetLepTrackEta, jetLepTrackPhi, pfMET, pfMETPhi;
	static Int_t nVtx;
	static Float_t *jetMt, *jetNHF, *jetCHF, *jetNEF, *jetVtxPt, *jetVtxMass, *jetVtx3DVal, *jetVtxNtrks, *jetVtx3DSig;
	static Int_t *jetNCH;
	
	jetPt                    = thisJet->Pt();
	jetEn                    = thisJet->Et();
	jetEta                   = thisJet->Eta();
	jetPhi                   = thisJet->Phi();
	jetMt                    = data.GetPtrFloat("THINjetMt");
	nVtx                     = data.GetInt("nVtx");
	jetRawPt                 = RawJet->Pt();
	jetLeadTrackPt           = leadJet->Pt();

	jetLepTrackPt            = lepJet->Pt();
	jetLepTrackEta           = lepJet->Eta();
	jetLepTrackPhi           = lepJet->Phi();

	jetNHF                   = data.GetPtrFloat("THINjetNHadEF");
	jetCHF                   = data.GetPtrFloat("THINjetCHadEF");
	jetNEF                   = data.GetPtrFloat("THINjetNEmEF");
	jetNCH                   = data.GetPtrInt("THINjetCMulti");
	jetVtxPt                 = data.GetPtrFloat("THINjetVtxPt");
	jetVtxMass               = data.GetPtrFloat("THINjetVtxMass");
	jetVtx3DVal              = data.GetPtrFloat("THINjetVtx3DVal");
	jetVtxNtrks              = data.GetPtrFloat("THINjetVtxNtrks");
	jetVtx3DSig              = data.GetPtrFloat("THINjetVtx3DSig");
	pfMET                    = data.GetFloat("pfMetCorrPt");
	pfMETPhi                 = data.GetFloat("pfMetCorrPhi");
  	// classification variables
	static float Jet_pt_, Jet_corr_,  Jet_eta_, Jet_mt_;
	static float Jet_leadTrackPt_, Jet_leptonPtRel_, Jet_leptonPt_, Jet_leptonDeltaR_;
	static float Jet_neHEF_, Jet_neEmEF_, Jet_chMult_, Jet_vtxPt_, Jet_vtxMass_, Jet_vtx3dL_, Jet_vtxNtrk_, Jet_vtx3deL_;
	static float Jet_PFMET_, Jet_METDPhi_, Jet_JetDR_;
	static float nVtx_;
  	// MVA classifiers for b-jet
	static TMVA::Reader* tmvaReader = NULL;
  	// one-time MVA initialization
	if (!tmvaReader) {
		tmvaReader = new TMVA::Reader("!Color:Silent");
		tmvaReader->AddVariable("Jet_pt",             	&Jet_pt_);
		tmvaReader->AddVariable("Jet_eta",            	&Jet_eta_);
		tmvaReader->AddVariable("Jet_mt",             	&Jet_mt_);
		tmvaReader->AddVariable("Jet_leadTrackPt",    	&Jet_leadTrackPt_);
		tmvaReader->AddVariable("Jet_leptonPtRel_new",  &Jet_leptonPtRel_);
		tmvaReader->AddVariable("Jet_leptonPt",       	&Jet_leptonPt_);
		tmvaReader->AddVariable("Jet_leptonDeltaR",   	&Jet_leptonDeltaR_);
		tmvaReader->AddVariable("Jet_neHEF",			      &Jet_neHEF_);
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
 	} // one-time initialization
  	// set MVA variables
  	Jet_pt_           = jetPt;//Jet pT
  	Jet_corr_         = jetPt/jetRawPt;//JEC
  	nVtx_             = nVtx;//nPVs
  	Jet_eta_          = jetEta;//Jet η
  	Jet_mt_           = jetMt[i];//Jet transverse mass
  	Jet_leadTrackPt_  = (jetLeadTrackPt < 0) ? 0. : jetLeadTrackPt;//pTLeadTrk, transverse momentum of the leading track in the jet
  	TLorentzVector TJet, TJetLep;
  	//TJet.SetPtEtaPhiE(jetPt,jetEta,jetPhi,jetEn);
  	TJetLep.SetPtEtaPhiE(jetLepTrackPt,jetLepTrackEta,jetLepTrackPhi,jetLepTrackPt*cosh(jetLepTrackEta));
  	//TVector3 TJetAxis(TJet.Vect().X(),TJet.Vect().Y(),TJet.Vect().Z());
  	TVector3 TJetAxis = thisJet->Vect();
  	Jet_leptonPtRel_	= (jetLepTrackPt > 0 && jetLepTrackPt<99900) ?  TJetLep.Perp(TJetAxis)	:0.;//Soft Lepton pTRel, relative transverse momentum of soft lepton candidate in the jet;
  	Jet_leptonPt_     = (jetLepTrackPt < 0 || jetLepTrackPt>99900) ? 0. : jetLepTrackPt;//Soft Lepton pT, transverse momentum of soft lepton candidate in the jet;
  	Jet_leptonDeltaR_ = (jetLepTrackPt < 0 || jetLepTrackPt>99900) ? 0. : deltaR(jetEta, jetPhi, jetLepTrackEta, jetLepTrackPhi);//Soft Lepton dR, relative distance in the η-phi space of soft lepton candidate in the jet and the jet;
  	Jet_neHEF_       	= (jetNHF[i]<1.) ? jetNHF[i]:1.;//Neutral hadron energy fraction
  	Jet_neEmEF_      	= (jetNEF[i]<1.) ? jetNEF[i]:1.;//Photon energy fraction
  	//Jet_chMult_       = (Float_t) jetNCH[i];//total number of jet constituents
  	Jet_vtxPt_        = jetVtxPt[i];//SecVtxPt, pT of the jet secondary vertex
  	Jet_vtxMass_      = jetVtxMass[i];//SecVtxM, Mass of the jet secondary vertex
  	Jet_vtx3dL_      	= (jetVtx3DVal[i]>0) ? jetVtx3DVal[i]:0.;//SecVtxdL, the 3-d flight length of the jet secondary vertex
  	Jet_vtxNtrk_      = (Int_t) jetVtxNtrks[i];//SecVtxNtrk, track multiplicity of the reconstructed secondary vertex
  	Jet_vtx3deL_		= (jetVtx3DSig[i] >0)		?  Jet_vtx3dL_/jetVtx3DSig[i]	:0.;//SecVtxdeL, Error on the 3-d flight length of the jet secondary vertex
  	Jet_PFMET_        = pfMET;
  	Jet_METDPhi_      = deltaPhi(pfMETPhi,jetPhi);
  	Jet_JetDR_        = jjDR;
  	return tmvaReader->EvaluateRegression("BDTG method")[0];
}


struct HA0JetInfo
{
	TLorentzVector* Hjet1   = NULL;
	TLorentzVector* Hjet2   = NULL;
	TLorentzVector* A0jet1  = NULL;
	TLorentzVector* A0jet2  = NULL;
	
	HA0JetInfo *weight4Jets = NULL;
private:
	float weight[4] = {-99,-99,-99,-99};
	int   ind[4]    = {-1,-1,-1,-1};
public:
	HA0JetInfo(){
		Hjet1       = NULL;
		Hjet2       = NULL;
		A0jet1      = NULL;
		A0jet2      = NULL;
		weight4Jets = NULL;
	}
	HA0JetInfo(TMVAinputInfo &info, int hi1, int hi2, int a0i1, int a0i2, bool doWeight = false) {
		ind[0] = hi1;
		ind[1] = hi2;
		ind[2] = a0i1;
		ind[3] = a0i2;
		Hjet1 = info.ak4jet[hi1];
		Hjet2 = info.ak4jet[hi2];
		A0jet1 = info.ak4jet[a0i1];
		A0jet2 = info.ak4jet[a0i2];
		if (doWeight) {
			if (weight[0]<0) weight[0] = TMVA_15plus3_jetGenJet_nu_3_1(info, hi1,Hjet1->DeltaR(*Hjet2));
			if (weight[1]<0) weight[1] = TMVA_15plus3_jetGenJet_nu_3_1(info, hi2,Hjet1->DeltaR(*Hjet2));
		}
	}
	~HA0JetInfo(){
		if (weight4Jets) {
			if (weight[0]>0) delete weight4Jets->Hjet1;
			if (weight[1]>0) delete weight4Jets->Hjet2;
			if (weight[2]>0) delete weight4Jets->A0jet1;
			if (weight[3]>0) delete weight4Jets->A0jet2;
			delete weight4Jets;
			weight4Jets = NULL;
		}
	}
	HA0JetInfo* weightJets(){
		if(!weight4Jets) weight4Jets = new HA0JetInfo();
		if (weight[0]<0||weight[1]<0) cout << "does not get weight" << endl;
		if (weight[0]>0) {
			weight4Jets->Hjet1 = new TLorentzVector();
			weight4Jets->Hjet1->SetPtEtaPhiE(Hjet1->Pt()*weight[0],Hjet1->Eta(),Hjet1->Phi(),Hjet1->E()*weight[0]);
		}
		else weight4Jets->Hjet1 = Hjet1;
		if (weight[1]>0) {
			weight4Jets->Hjet2 = new TLorentzVector();
			weight4Jets->Hjet2->SetPtEtaPhiE(Hjet2->Pt()*weight[1],Hjet2->Eta(),Hjet2->Phi(),Hjet2->E()*weight[1]);
		}
		else weight4Jets->Hjet2 = Hjet2;
		weight4Jets->A0jet1 = A0jet1;
		weight4Jets->A0jet2 = A0jet2;
		weight4Jets->weight[0] = -1;
		weight4Jets->weight[1] = -1;
		weight4Jets->weight[2] = -1;
		weight4Jets->weight[3] = -1;
		weight4Jets->ind[0] = ind[0];
		weight4Jets->ind[1] = ind[1];
		weight4Jets->ind[2] = ind[2];
		weight4Jets->ind[3] = ind[3];
		// cout << weight4Jets->weight[0] << "\t" << weight4Jets->weight[1] << "\t" << weight4Jets->weight[2] << "\t" << weight4Jets->weight[3] << "\t" << endl;
		return weight4Jets;
	}
	inline float Mh()     	{return (*Hjet2+*Hjet1).M();}
	inline float MA0()    	{return (*A0jet2+*A0jet1).M();}
	inline float MZp()    	{return (*Hjet2+*Hjet1+*A0jet2+*A0jet1).M();}
	inline float Pth()	  	{return (*Hjet2+*Hjet1).Pt();}
	inline float Pta0()		{return (*A0jet2+*A0jet1).Pt();}
	inline float PtZp()		{return (*Hjet2+*Hjet1+*A0jet2+*A0jet1).Pt();}
	inline float HdR()    	{return Hjet1->DeltaR(*Hjet2);}
	inline float A0dR()   	{return A0jet1->DeltaR(*A0jet2);}
	inline float ZpdR()   	{return (*Hjet2+*Hjet1).DeltaR(*A0jet2+*A0jet1);}
	inline float xSqure() 	{return calChiSquare(Mh(),MA0(),MZp());}
	inline float* getWeight()      	{return weight;}
	inline float  getWeight(int i) 	{return (i<4&&i>=0)?weight[i]:-1;}
	inline int* getIndex()			{return ind;}
	inline int  getIndex(int i)  	{return (i<4&&i>=0)?ind[i]:-1;}
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
};
