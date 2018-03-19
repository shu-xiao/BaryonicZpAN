//#include "Utilities.h"
#include "TMVA/Reader.h"


//2017_3_1 update==============
float deltaR(float eta1,float phi1, float eta2,float phi2) {
	return sqrt( (eta2-eta1)*(eta2-eta1) + (phi1-phi2)*(phi1-phi2) );
}
float deltaPhi(float phi1, float phi2) {
	return (abs(phi2-phi1)<(Float_t)TMath::Pi())? abs(phi2-phi1): 2*TMath::Pi()-abs(phi2-phi1);
}
Float_t TMVA_15plus3_jetGenJet_nu_3_1(TreeReader &data, Int_t i, Float_t jjDR) 
{
	static TClonesArray* ak4jetP4, *rawjetP4, *THINjetLeadTrackP4, *LepjetP4;
	static Int_t nJet;
	static TLorentzVector* RawJet, *thisJet, *leadJet, *lepJet;
	ak4jetP4            =  (TClonesArray*) data.GetPtrTObject("THINunCorrJetP4");
	rawjetP4            =  (TClonesArray*) data.GetPtrTObject("THINjetP4");
	THINjetLeadTrackP4  =  (TClonesArray*) data.GetPtrTObject("THINjetLeadTrackP4");
	LepjetP4            =  (TClonesArray*) data.GetPtrTObject("THINjetLepTrackP4");
	nJet                     = data.GetInt("THINnJet");
	
	Float_t jetPtList[nJet], jetEnList[nJet] ;
	Float_t jetEtaList[nJet], jetPhiList[nJet];
	Float_t jetRawPtList[nJet];
	Float_t jetLeadTrackPtList[nJet];
	Float_t jetLepTrackPtList[nJet], jetLepTrackEtaList[nJet], jetLepTrackPhiList[nJet];
	for (int i=0;i<nJet;i++) {
		RawJet  = (TLorentzVector*)rawjetP4->At(i);
		jetRawPtList[i] = RawJet->Pt(); 
		
		thisJet = (TLorentzVector*)ak4jetP4->At(i);
		jetPtList[i]  = thisJet->Pt();
		jetEnList[i]  = thisJet->Et();		
		jetEtaList[i] = thisJet->Eta();
		jetPhiList[i] = thisJet->Phi();
		
		leadJet = (TLorentzVector*)THINjetLeadTrackP4->At(i);
		jetLeadTrackPtList[i] = leadJet->Pt();
		
		lepJet = (TLorentzVector*)LepjetP4->At(i);
		jetLepTrackPtList[i]  = lepJet->Pt();
		jetLepTrackPhiList[i] = lepJet->Phi();
		jetLepTrackEtaList[i] = lepJet->Eta();

	}
  	//Int_t    nVtx                     = data.GetInt("nGoodVtx");
  	//Float_t* jetPt                    = data.GetPtrFloat("jetPt");
    //Float_t* jetEn                    = data.GetPtrFloat("jetEn");
    //Float_t* jetEta                   = data.GetPtrFloat("jetEta");
    //Float_t* jetPhi                   = data.GetPtrFloat("jetPhi");
    //Float_t* jetMt                    = data.GetPtrFloat("jetMt");
	Float_t* jetPt                    = jetPtList;
	Float_t* jetEn                    = jetEnList;
	Float_t* jetEta                   = jetEtaList;
	Float_t* jetPhi                   = jetPhiList;
	Float_t* jetMt                    = data.GetPtrFloat("THINjetMt");
	Int_t    nVtx                     = data.GetInt("nVtx");
	// Float_t* jetRawPt                 = data.GetPtrFloat("jetRawPt");
	Float_t* jetRawPt                 = jetRawPtList;
	Float_t* jetLeadTrackPt           = jetLeadTrackPtList;
	// Float_t* jetLeadTrackPt           = data.GetPtrFloat("jetLeadTrackPt");
	// Float_t* jetLepTrackEta           = data.GetPtrFloat("jetLepTrackEta");
	// Float_t* jetLepTrackPhi           = data.GetPtrFloat("jetLepTrackPhi");
	// Float_t* jetLepTrackPt            = data.GetPtrFloat("jetLepTrackPt");
	Float_t* jetLepTrackPt            = jetLepTrackPtList;
	Float_t* jetLepTrackEta           = jetLepTrackEtaList;
	Float_t* jetLepTrackPhi           = jetLepTrackPhiList;
	// Float_t* jetNHF                   = data.GetPtrFloat("jetNHF");
	// Float_t* jetCHF                   = data.GetPtrFloat("jetCHF");
	// Float_t* jetNEF                   = data.GetPtrFloat("jetNEF");
	// Int_t*   jetNCH                   = data.GetPtrInt("jetNCH");
	// Float_t* jetVtxPt                 = data.GetPtrFloat("jetVtxPt");
	// Float_t* jetVtxMass               = data.GetPtrFloat("jetVtxMass");
	// Float_t* jetVtx3DVal              = data.GetPtrFloat("jetVtx3DVal");
	// Float_t* jetVtxNtrks              = data.GetPtrFloat("jetVtxNtrks");
	// Float_t* jetVtx3DSig              = data.GetPtrFloat("jetVtx3DSig");
	// Float_t  pfMET                    = data.GetFloat("pfMET");
	// Float_t  pfMETPhi                 = data.GetFloat("pfMETPhi");
	Float_t* jetNHF                   = data.GetPtrFloat("THINjetNHadEF");
	Float_t* jetCHF                   = data.GetPtrFloat("THINjetCHadEF");
	Float_t* jetNEF                   = data.GetPtrFloat("THINjetNEmEF");
	Int_t*   jetNCH                   = data.GetPtrInt("THINjetCMulti");
	Float_t* jetVtxPt                 = data.GetPtrFloat("THINjetVtxPt");
	Float_t* jetVtxMass               = data.GetPtrFloat("THINjetVtxMass");
	Float_t* jetVtx3DVal              = data.GetPtrFloat("THINjetVtx3DVal");
	Float_t* jetVtxNtrks              = data.GetPtrFloat("THINjetVtxNtrks");
	Float_t* jetVtx3DSig              = data.GetPtrFloat("THINjetVtx3DSig");
	Float_t  pfMET                    = data.GetFloat("pfMetCorrPt");
	Float_t  pfMETPhi                 = data.GetFloat("pfMetCorrPhi");
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
  	Jet_pt_           = jetPt[i];//Jet pT
  	Jet_corr_         = jetPt[i]/jetRawPt[i];//JEC
  	nVtx_             = nVtx;//nPVs
  	Jet_eta_          = jetEta[i];//Jet η
  	Jet_mt_           = jetMt[i];//Jet transverse mass
  	Jet_leadTrackPt_  = (jetLeadTrackPt[i] < 0) ? 0. : jetLeadTrackPt[i];//pTLeadTrk, transverse momentum of the leading track in the jet
  	TLorentzVector TJet, TJetLep;
  	TJet.SetPtEtaPhiE(jetPt[i],jetEta[i],jetPhi[i],jetEn[i]);
  	TJetLep.SetPtEtaPhiE(jetLepTrackPt[i],jetLepTrackEta[i],jetLepTrackPhi[i],jetLepTrackPt[i]*cosh(jetLepTrackEta[i]));
  	TVector3 TJetAxis(TJet.Vect().X(),TJet.Vect().Y(),TJet.Vect().Z());
  	Jet_leptonPtRel_	= (jetLepTrackPt[i] > 0) ?  TJetLep.Perp(TJetAxis)	:0.;//Soft Lepton pTRel, relative transverse momentum of soft lepton candidate in the jet;
  	Jet_leptonPt_     = (jetLepTrackPt[i] < 0) ? 0. : jetLepTrackPt[i];//Soft Lepton pT, transverse momentum of soft lepton candidate in the jet;
  	Jet_leptonDeltaR_ = (jetLepTrackPt[i] < 0) ? 0. : deltaR(jetEta[i], jetPhi[i], jetLepTrackEta[i], jetLepTrackPhi[i]);//Soft Lepton dR, relative distance in the η-phi space of soft lepton candidate in the jet and the jet;
  	Jet_neHEF_       	= (jetNHF[i]<1.) ? jetNHF[i]:1.;//Neutral hadron energy fraction
  	Jet_neEmEF_      	= (jetNEF[i]<1.) ? jetNEF[i]:1.;//Photon energy fraction
  	//Jet_chMult_       = (Float_t) jetNCH[i];//total number of jet constituents
  	Jet_vtxPt_        = jetVtxPt[i];//SecVtxPt, pT of the jet secondary vertex
  	Jet_vtxMass_      = jetVtxMass[i];//SecVtxM, Mass of the jet secondary vertex
  	Jet_vtx3dL_      	= (jetVtx3DVal[i]>0) ? jetVtx3DVal[i]:0.;//SecVtxdL, the 3-d flight length of the jet secondary vertex
  	Jet_vtxNtrk_      = (Int_t) jetVtxNtrks[i];//SecVtxNtrk, track multiplicity of the reconstructed secondary vertex
  	Jet_vtx3deL_		= (jetVtx3DSig[i] >0)		?  Jet_vtx3dL_/jetVtx3DSig[i]	:0.;//SecVtxdeL, Error on the 3-d flight length of the jet secondary vertex
  	Jet_PFMET_        = pfMET;
  	Jet_METDPhi_      = deltaPhi(pfMETPhi,jetPhi[i]);
  	Jet_JetDR_        = jjDR;
  	return tmvaReader->EvaluateRegression("BDTG method")[0];
}

