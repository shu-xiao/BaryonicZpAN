#include "pfWeight.h"



using namespace std;

const float CISVV2_CUT_L = 0.5426;
const float CISVV2_CUT_M = 0.8484;
const float CISVV2_CUT_T = 0.9535;
const int maxCom = 50;
Float_t CISVV2, CISVV2_1[maxCom], CISVV2_2[maxCom];
Float_t hPt[maxCom],Mh[maxCom],hDeltaR[maxCom],hDeltaEta[maxCom], hDeltaPhi[maxCom],hptas[maxCom],hsdas[maxCom];
bool defaultSel(int i=0,const float& btag=CISVV2_CUT_L) {return true;}
bool minSel(int i, const float& btag=CISVV2_CUT_L) {return min(CISVV2_1[i],CISVV2_2[i])>=btag;}
bool maxSel(int i, const float& btag=CISVV2_CUT_L) {return max(CISVV2_1[i],CISVV2_2[i])>=btag;}
bool highPtSel(int i, const float& btag=CISVV2_CUT_L) {return CISVV2_1[i]>=btag;}
bool lowPtSel(int i,const float& btag=CISVV2_CUT_L) {return CISVV2_2[i]>=btag;}
bool meanSel(int i,const float& btag=CISVV2_CUT_L) {return (CISVV2_1[i]+CISVV2_2[i])/2>=btag;}




void pfMain(bool dotest=true){
	//pfWeight a;
	//a.openFile();
	//ClassImp(nHist);
	//cout << "start" << endl;
	dotest = 0;
	const string fName[9] = {"tree_QCD_TMVA_HT50to100.root","tree_QCD_TMVA_HT100to200.root","tree_QCD_TMVA_HT200to300.root","tree_QCD_TMVA_HT300to500.root","tree_QCD_TMVA_HT500to700.root","tree_QCD_TMVA_HT700to1000.root","tree_QCD_TMVA_HT1000to1500.root","tree_QCD_TMVA_HT1500to2000.root","tree_QCD_TMVA_HT2000toInf.root"};

	const float xsHTbeam[9] = {246400000,27990000,1712000,347700,32100,6831,1207,119.9,25.24};
    const float L2016=35.9*1000;//35.9 fb^-1
    const vector<string> selstr = {"MinBTag","MaxBTag","MeanBTag","HighPt","LowPt"};
    const int  nSel = 5;
    Int_t nCom, nEvents = 1;
    TFile *f = NULL, *fWrite;
    fWrite = new TFile("diffvar.root","recreate");
    TTree *t1 = NULL;
    TCanvas *c1 = new TCanvas("c1","c1",3);
    //bool (*sel_1)(float ,const float&)=0;
    //bool (*sel_2)(float ,const float&)=0;

	nHist diffSel[nSel];

	for (int i=0;i<nSel;i++) diffSel[i] = nHist(selstr[i]);
	//cout << "ini.." << endl;
	bool (*nSelFun[nSel])(int i,const float& btag)={minSel,maxSel,meanSel,highPtSel,lowPtSel};
	float b = 1;
	vector<vector<int>> failIndex(9*nSel);
	// loop for different HT beam samples
	for (int xs=0;xs<9;xs++) {
		
		if (xs==0) continue;
		cout << "Processing " << fName[xs] << endl;
		f = new TFile(fName[xs].data());
		TH1F* h_allEvent = (TH1F*)f->Get("h_allEvent");
        nEvents = h_allEvent->GetEntries();
		if (!dotest) b = L2016*xsHTbeam[xs]/nEvents;


		t1 = (TTree*)f->Get("tree");
		t1->SetBranchAddress("nCom",&nCom);
        t1->SetBranchAddress("CISVV2_Hb1",CISVV2_1);
        t1->SetBranchAddress("CISVV2_Hb2",CISVV2_2);
        t1->SetBranchAddress("Mh_weightLT",Mh);
        t1->SetBranchAddress("hPt",hPt);
        t1->SetBranchAddress("hDeltaR",hDeltaR);
        t1->SetBranchAddress("hDeltaPhi",hDeltaPhi);
        t1->SetBranchAddress("hDeltaEta",hDeltaEta);
        t1->SetBranchAddress("hptAs",hptas);
        t1->SetBranchAddress("hsdAs",hsdas);

		static int i,j, k;
		for (i=0;i<nSel;i++) failIndex[xs*nSel+i] = std::vector<int> (t1->GetEntries(),-1);
		// loop for events
		// get pass and fail events
		for (i=0;i<t1->GetEntries();i++) {
			t1->GetEntry(i);
			bool isPass[nSel] = {0};
        	int ind[nSel] = {0};
        	int find[nSel];
        	for (int kk=0;kk<nSel;kk++) find[kk] = -1;
        	// combination loop for pass events
			for (j=0;j<nCom;j++) {
				if (Mh[j]<90||Mh[j]>160) continue;
				// loop for different var
				for (k=0;k<nSel;k++) {
					if (!isPass[k]&&nSelFun[k](j,CISVV2_CUT_L)){
						isPass[k] = true;
						ind[k] = j;
						failIndex[xs*nSel+k][i]=99999;
					}
				}
			}
			// combination loop for fail event
			for (k=0;k<nSel;k++) { 
				for (j=0;j<nCom;j++) {
					if (Mh[j]<90||Mh[j]>160) continue;
					if (!isPass[k]&&!nSelFun[k](j,CISVV2_CUT_L)) {
						find[k] = j;
						failIndex[xs*nSel+k][i]=j;
						break;
					}
				}
			}
			// fill hist with different var
			for (k=0;k<nSel;k++) {
				if (!isPass[k]&&find[k]<0) continue;
				if (find[k]>=0&&!isPass[k]) ind[k] = find[k];
				diffSel[k].h_Mh[i%2*2+isPass[k]]->Fill(Mh[ind[k]],b);
				diffSel[k].h_hPt[i%2*2+isPass[k]]->Fill(hPt[ind[k]],b);
				diffSel[k].h_hDeltaR[i%2*2+isPass[k]]->Fill(hDeltaR[ind[k]],b);
				diffSel[k].h_hDeltaEta[i%2*2+isPass[k]]->Fill(hDeltaEta[ind[k]],b);
				diffSel[k].h_hDeltaPhi[i%2*2+isPass[k]]->Fill(hDeltaPhi[ind[k]],b);
				diffSel[k].h_hsdas[i%2*2+isPass[k]]->Fill(hsdas[ind[k]],b);
				diffSel[k].h_hptas[i%2*2+isPass[k]]->Fill(hptas[ind[k]],b);
			}
		} // end of loop events

		f->Close();
		f = NULL;
		if (dotest) break;
	}// end of loop HT beam

	cout << "cal pf ratio" << endl;
	for (int i=0;i<nSel;i++) {
		diffSel[i].divide();
		diffSel[i].fit();
	}
	doReject=false;
	// loop for different HT beam samples
	// weight
	for (int xs=0;xs<9;xs++) {
		if (xs==0) continue;
		cout << "Processing " << fName[xs] << endl;
		f = new TFile(fName[xs].data());
		TH1F* h_allEvent = (TH1F*)f->Get("h_allEvent");
		
        nEvents = h_allEvent->GetEntries();
		if (!dotest) b = L2016*xsHTbeam[xs]/nEvents;
		t1 = (TTree*)f->Get("tree");
		t1->SetBranchAddress("nCom",&nCom);
        t1->SetBranchAddress("CISVV2_Hb1",CISVV2_1);
        t1->SetBranchAddress("CISVV2_Hb2",CISVV2_2);
        t1->SetBranchAddress("Mh_weightLT",Mh);
        t1->SetBranchAddress("hPt",hPt);
        t1->SetBranchAddress("hDeltaR",hDeltaR);
        t1->SetBranchAddress("hDeltaPhi",hDeltaPhi);
        t1->SetBranchAddress("hDeltaEta",hDeltaEta);
        t1->SetBranchAddress("hptAs",hptas);
        t1->SetBranchAddress("hsdAs",hsdas);
		
		static int i,j,k;
		static float bw[2];
		// fill hist with different var
		for (i=0;i<t1->GetEntries();i++) {
			if (i%2==0) continue;
			t1->GetEntry(i);
			bool isPass[nSel] = {0};
        	int ind[nSel] = {0};
        	int find[nSel];
        	for (int kk=0;kk<nSel;kk++) find[kk] = -1;
            for (int j=0;j<nCom;j++) {
                
                if (Mh[j]<90||Mh[j]>160) continue;
                for (k=0;k<nSel;k++) {
                	if (!isPass[k]&&nSelFun[k](j,CISVV2_CUT_L)) {
                    	isPass[k] = true;
                    	ind[k] = j;
                	}
            	}
            }
            // veto pass, leave fail
            for (k=0;k<nSel;k++) { 
            	for (int j=0;j<nCom;j++) {
                	if (Mh[j]<90||Mh[j]>160) continue;
                	if (!isPass[k]&&!nSelFun[k](j,CISVV2_CUT_L)) {
                    	find[k] = j;
                    	break;
                	}
                }
            }
            //if (find<0) continue;

			// loop for selection
			for (k=0;k<nSel;k++) {
				if (find[k]<0) continue;
				// loop for different fit function
				//find = failIndex[xs*nSel+k][i];
				bw[0] = b*diffSel[k].fitfun1->Eval(Mh[find[k]]);
				bw[1] = b*diffSel[k].fitfun2->Eval(Mh[find[k]]);
				
				for (j=0;j<2;j++) {
        			diffSel[k].h_Mh[5+j]->Fill(Mh[find[k]],bw[j]);
           			diffSel[k].h_hPt[5+j]->Fill(hPt[find[k]],bw[j]);
            		diffSel[k].h_hDeltaR[5+j]->Fill(hDeltaR[find[k]],bw[j]);
            		diffSel[k].h_hDeltaPhi[5+j]->Fill(hDeltaPhi[find[k]],bw[j]);
            		diffSel[k].h_hDeltaEta[5+j]->Fill(hDeltaEta[find[k]],bw[j]);
            		diffSel[k].h_hptas[5+j]->Fill(hptas[find[k]],bw[j]);
            		diffSel[k].h_hsdas[5+j]->Fill(hsdas[find[k]],bw[j]);
            	}
            }
			// error code
			/*
			bool isPass[nSel] = {0};
			// loop for selection
			for (k=0;k<nSel;k++) {
				if (failIndex[xs*nSel+k][i]>=maxCom||failIndex[xs*nSel+k][i]<0) isPass[k] = true;
				if (isPass[k]) continue;
				
				// loop for different fit function
				find = failIndex[xs*nSel+k][i];
				bw[0] = b*diffSel[k].fitfun1->Eval(Mh[find]);
				bw[1] = b*diffSel[k].fitfun2->Eval(Mh[find]);
				
				for (j=0;j<2;j++) {
        			diffSel[k].h_Mh[5+j]->Fill(Mh[find],bw[j]);
           			diffSel[k].h_hPt[5+j]->Fill(hPt[find],bw[j]);
            		diffSel[k].h_hDeltaR[5+j]->Fill(hDeltaR[find],bw[j]);
            		diffSel[k].h_hDeltaPhi[5+j]->Fill(hDeltaPhi[find],bw[j]);
            		diffSel[k].h_hDeltaEta[5+j]->Fill(hDeltaEta[find],bw[j]);
            		diffSel[k].h_hptas[5+j]->Fill(hptas[find],bw[j]);
            		diffSel[k].h_hsdas[5+j]->Fill(hsdas[find],bw[j]);
            	}
        	}
        	*/

		}// end of loop of entries


		f->Close();
		f = NULL;
		if (dotest) break;
	}// end of loop HT beam
	cout << "end of weight" << endl;
	//for (int i=0;i<nSel;i++) diffSel[i].drawDiffAll();
	//string nTitle[7] = {"M_{h}","higgs Pt (GeV)","#DeltaR_{bb}","#Delta#eta_{bb}","#Delta#phi_{bb}",""};
	
	drawDiff(diffSel[0].h_Mh[5],diffSel[0].h_Mh[3],"M_{h}",-1);
	for (int i=0;i<nSel;i++)  for (int j=5;j<=6;j++) drawDiff(diffSel[i].h_Mh[j],diffSel[i].h_Mh[3],"M_{h}");
	for (int i=0;i<nSel;i++)  for (int j=5;j<=6;j++) drawDiff(diffSel[i].h_hPt[j],diffSel[i].h_hPt[3],"higgs Pt (GeV)");
	for (int i=0;i<nSel;i++)  for (int j=5;j<=6;j++) drawDiff(diffSel[i].h_hDeltaR[j],diffSel[i].h_hDeltaR[3],"#DeltaR_{bb}");
	for (int i=0;i<nSel;i++)  for (int j=5;j<=6;j++) drawDiff(diffSel[i].h_hDeltaEta[j],diffSel[i].h_hDeltaEta[3],"#Delta#eta_{bb}");
	for (int i=0;i<nSel;i++)  for (int j=5;j<=6;j++) drawDiff(diffSel[i].h_hDeltaPhi[j],diffSel[i].h_hDeltaPhi[3],"#Delta#phi_{bb}");
	for (int i=0;i<nSel;i++)  for (int j=5;j<=6;j++) drawDiff(diffSel[i].h_hptas[j],diffSel[i].h_hptas[3],"pt assymetry (min(Pt1,Pt2)#Delta R/m_{jj})^{2}");
	for (int i=0;i<nSel;i++)  for (int j=5;j<=6;j++) drawDiff(diffSel[i].h_hsdas[j],diffSel[i].h_hsdas[3],"min(pt1,pt2)/(pt1+pt2)");
	
	drawDiff(diffSel[0].h_Mh[5],diffSel[0].h_Mh[3],"M_{h}",-2);
	
	//drawDiff(diffSel[0].h_Mh[5],diffSel[0].h_Mh[3],"M_{h}",1);
    //drawDiff(diffSel[0].h_Mh[6],diffSel[0].h_Mh[3],"M_{h}",2);
	fWrite->Write();
	fWrite->Close();
}
