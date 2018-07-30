

using namespace std;
bool doReject=true;
int nVar = 0;
bool isDraw = false;
string pdfName="pfRatioCompareHead.pdf";
void setMax(TH1F* h1, TH1F* h2) {
    float max = h1->GetMaximum();
    float max2 = h2->GetMaximum();
    if (max<max2) max = max2;
    h1->SetMaximum(max*1.3);
    h2->SetMaximum(max*1.3);
}
void histSetting(TH1F* hist,int i=0){
	hist->Sumw2();
	hist->SetLineWidth(2);
	if (i==5) hist->SetLineColor(kRed);
	else if (i==6) hist->SetLineColor(kCyan);
}
double linear(double *x,double *par) {
    if (doReject&&x[0]<140&&x[0]>120) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0]+par[1]*x[0];
}
double pol_2(double *x,double *par) {
    if (doReject&&x[0]<140&&x[0]>120) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}


const vector<string> suf = {"_fail","_pass","_testF","_testP","_PFratio","_weight","_weight2"};

struct nHist {
	
	static const int nHistlen = 7;
	/*
	static int nVar;
	static string pdfName;
	static bool isDraw;
	*/
	bool isIni = false;
	TF1*  fitfun1,*fitfun2;

    TH1F* h_Mh[nHistlen]={0};
    TH1F* h_hPt[nHistlen]={0};
    TH1F* h_hDeltaR[nHistlen]={0};
    TH1F* h_hDeltaEta[nHistlen]={0};
    TH1F* h_hDeltaPhi[nHistlen]={0};
    TH1F* h_hptas[nHistlen]={0}, *h_hsdas[nHistlen]={0};
    TCanvas *c3;
    nHist(){}
    nHist(string sel) {
    	//isIni = true;
    	//cout << "con." << endl;
    	//cout << "con...." << endl;
	    for (int i=0;i<nHistlen;i++) {
	        h_Mh[i] = new TH1F(Form("h_Mh_%s%s",sel.data(),suf[i].data()),Form("h_Mh_%s%s",sel.data(),suf[i].data()),20,90,160);
	        h_hPt[i] = new TH1F(Form("h_hPt_%s%s",sel.data(),suf[i].data()),Form("h_hPt_%s%s",sel.data(),suf[i].data()),60,0,900);
	        h_hDeltaR[i] = new TH1F(Form("h_hDeltaR_%s%s",sel.data(),suf[i].data()),Form("h_hDeltaR_%s%s",sel.data(),suf[i].data()),20,0,4);
	        h_hDeltaPhi[i] = new TH1F(Form("h_hDeltaPhi_%s%s",sel.data(),suf[i].data()),Form("h_hDeltaPhi_%s%s",sel.data(),suf[i].data()),32,0,3.2);
	        h_hDeltaEta[i] = new TH1F(Form("h_hDeltaEta_%s%s",sel.data(),suf[i].data()),Form("h_hDeltaEta_%s%s",sel.data(),suf[i].data()),20,0,4);
	        h_hptas[i] = new TH1F(Form("h_hptas_%s%s",sel.data(),suf[i].data()),Form("h_hptas_%s%s",sel.data(),suf[i].data()),20,0,2);
	        h_hsdas[i] = new TH1F(Form("h_hsdas_%s%s",sel.data(),suf[i].data()),Form("h_hsdas_%s%s",sel.data(),suf[i].data()),30,0,0.6);
	        histSetting(h_Mh[i],i);
	        histSetting(h_hPt[i],i);
	        histSetting(h_hDeltaR[i],i);
	        histSetting(h_hDeltaPhi[i],i);
	        histSetting(h_hDeltaEta[i],i);
	        histSetting(h_hptas[i],i);
	        histSetting(h_hsdas[i],i);
	    }
	    fitfun1 = new TF1(Form("f1_%s",sel.data()),linear,90,160,2);
    	fitfun2 = new TF1(Form("f2_%s",sel.data()),pol_2,90,160,3);
	    fitfun1->SetLineWidth(2);
    	fitfun2->SetLineWidth(2);
    	fitfun1->SetLineColor(kBlue);
    	fitfun2->SetLineColor(kRed);
    	//cout << "con.........." << endl;
	}
	~nHist(){
		//if (nVar==0&&c3) c3->Print((pdfName+"]").data());
	}
	void draw() {}
	void release(){
		for (int i=0;i<nHistlen;i++) {
			delete h_Mh[i];			h_Mh[i]=0;
			delete h_hPt[i];		h_hPt[i]=0;
			delete h_hDeltaR[i];	h_hDeltaR[i]=0;
			delete h_hDeltaEta[i];	h_hDeltaEta[i]=0;
			delete h_hDeltaPhi[i];	h_hDeltaPhi[i]=0;
			delete h_hsdas[i];		h_hsdas[i]=0;
			delete h_hptas[i];		h_hptas[i]=0;
		}
	}
	void divide() {
		h_Mh[4]->Divide(h_Mh[1],h_Mh[0]);
    	h_hPt[4]->Divide(h_hPt[1],h_hPt[0]);
    	h_hDeltaR[4]->Divide(h_hDeltaR[1],h_hDeltaR[0]);
    	h_hDeltaEta[4]->Divide(h_hDeltaEta[1],h_hDeltaEta[0]);
    	h_hDeltaPhi[4]->Divide(h_hDeltaPhi[1],h_hDeltaPhi[0]);
    	h_hptas[4]->Divide(h_hptas[1],h_hptas[0]);
    	h_hsdas[4]->Divide(h_hsdas[1],h_hsdas[0]);
	}
	void fit() {
		h_Mh[4]->Fit(fitfun1);
    	h_Mh[4]->Fit(fitfun2,"+");
	}


};
/*
	void drawDiff(TH1F* h1, TH1F* h2, string title="") {
	    // h1 est, h2 MC
	    isDraw = true;
	    setMax(h1,h2);
	    if (!c3) {
	    	c3 = new TCanvas("c3","c3",3);
	    	c3->Divide(1,2,0.01,0.01);
	    	c3->GetPad(1)->SetLeftMargin(0.12);
	    	c3->GetPad(2)->SetLeftMargin(0.12);
	    	c3->GetPad(1)->SetBottomMargin(0.12);
	    	c3->GetPad(2)->SetGridy();
	    	c3->GetPad(2)->SetPad(0.0,0.0,1,0.3);
	    	c3->GetPad(1)->SetPad(0.0,0.3,1,1);
	    	c3->GetPad(1)->SetTicks();
	    	c3->GetPad(2)->SetTicks();
	    	c3->Print((pdfName+"[").data());
	    }
	    gStyle->SetOptStat(0);
	    gStyle->SetOptTitle(0);
	    c3->cd(1);
	    TH1F* h_copy = new TH1F(*h1);
	    h_copy->Divide(h2);
	    h_copy->SetLineWidth(2);
	    h1->SetLineColor(kBlue);
	    h2->SetLineColor(kBlack);
	    h1->GetYaxis()->SetTitle("A.U.");
	    h1->GetXaxis()->SetTitle(title.data());

	    h1->GetYaxis()->SetTitleSize(0.05);
	    h1->GetXaxis()->SetTitleSize(0.05);
	    //h1->GetXaxis()->SetTitleOffset(0.4);
	    //h1->GetYaxis()->SetTitleOffset(0.4);
	    h2->SetStats(0);
	    h1->Draw("e");
	    h2->Draw("esame");
	    // legend
	    TLegend leg(0.5,0.7,0.88,0.85);
	    TString htitle = h1->GetTitle();
	    //if (htitle.Contains("weight2")) leg.AddEntry(h1,"p/f(parabola funciton fit)*f");
	    //else leg.AddEntry(h1,"p/f(linear function fit)*f");
	    // set estimate title
	    if (htitle.Contains("weight2")) leg.SetHeader("Linear");
	    else leg.SetHeader("Parabola");
	    leg.AddEntry(h1,"Estimation");
	    leg.AddEntry(h2,"MC");
	    leg.SetBorderSize(0);
	    leg.Draw();
	    // ratio
	    c3->cd(2);

	    h_copy->GetYaxis()->SetRangeUser(0,2);
	    h_copy->GetYaxis()->SetTitle("Estimation/MC");
	    h_copy->GetYaxis()->CenterTitle();
	    h_copy->SetLineColor(kBlack);
	    h_copy->GetXaxis()->SetLabelSize(0);
	    h_copy->GetYaxis()->SetLabelSize(0.1);
	    h_copy->GetYaxis()->SetTitleSize(0.1);
	    h_copy->GetYaxis()->SetTitleOffset(0.5);
	    //h_copy->SetMarkerStyle(20);
	    // Y = 1 LINE
	    static TF1 *f1 = 0;
	    if (!f1) f1 = new TF1("f1","1",-1000,1000);
	    //f1->SetLineWidth(2);
	    h_copy->Draw("E1");
	    f1->Draw("same");
	    //h_copy->Draw("e1same");
	    
	    c3->Print(pdfName.data());
	    delete c3;
	    delete h_copy;
	}
	void drawDiffAll() {
		drawDiff(h_Mh[5],h_Mh[3],"M_{h}");
	    drawDiff(h_Mh[6],h_Mh[3],"M_{h}");
	    drawDiff(h_hPt[5],h_hPt[3],"higgs Pt");
	    drawDiff(h_hPt[6],h_hPt[3],"higgs Pt");
	    drawDiff(h_hDeltaR[5],h_hDeltaR[3],"#DeltaR_{bb}");
	    drawDiff(h_hDeltaR[6],h_hDeltaR[3],"#DeltaR_{bb}");
	    drawDiff(h_hDeltaEta[5],h_hDeltaEta[3], "#Delta#eta_{bb}");
	    drawDiff(h_hDeltaEta[6],h_hDeltaEta[3],"#Delta#eta_{bb}");
	    drawDiff(h_hDeltaPhi[5],h_hDeltaPhi[3],"#Delta#phi_{bb}");
	    drawDiff(h_hDeltaPhi[6],h_hDeltaPhi[3],"#Delta#phi_{bb}");
	    drawDiff(h_hptas[5],h_hptas[3],"pt assymetry (min(Pt1,Pt2)#Delta R/m_{jj})^{2}");
	    drawDiff(h_hptas[6],h_hptas[3],"pt assymetry (min(Pt1,Pt2)#Delta R/m_{jj})^{2}");
	    drawDiff(h_hsdas[5],h_hsdas[3],"min(pt1,pt2)/(pt1+pt2)");
	    drawDiff(h_hsdas[6],h_hsdas[3],"min(pt1,pt2)/(pt1+pt2)");
	}
	*/
/*
int nHist::nVar = 0;
bool nHist::isDraw = false;
string nHist::pdfName="pfRatioCompareHead.pdf";
*/

void drawDiff(TH1F* h1, TH1F* h2, string title="",int io=0, string fileName="pfRatioComparenVar.pdf") {
    // h1 est, h2 sim
    TCanvas* c4 = new TCanvas("c4","c4",3);
    if (io==-1) c4->Print((fileName+"[").data());
    else if (io==-2) c4->Print((fileName+"]").data());
    if (io<0) return;
    setMax(h1,h2);
    
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    c4->Divide(1,2,0.01,0.01);
    c4->cd(1);
    TH1F* h_copy = new TH1F(*h1);
    h_copy->Divide(h2);
    h1->SetLineWidth(2);
    h2->SetLineWidth(2);
    h_copy->SetLineWidth(2);
    h1->SetLineColor(kBlue);
    h2->SetLineColor(kBlack);
    h1->GetYaxis()->SetTitle("A.U.");
    h1->GetXaxis()->SetTitle(title.data());
    c4->GetPad(1)->SetLeftMargin(0.12);
    c4->GetPad(2)->SetLeftMargin(0.12);
    c4->GetPad(1)->SetBottomMargin(0.12);
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->GetXaxis()->SetTitleSize(0.05);
    //h1->GetXaxis()->SetTitleOffset(0.4);
    //h1->GetYaxis()->SetTitleOffset(0.4);
    h2->SetStats(0);
    h1->Draw("e");
    h2->Draw("esame");
    // legend
    TLegend leg(0.6,0.7,0.88,0.85);
    TString htitle = h1->GetTitle();
    //if (htitle.Contains("weight2")) leg.AddEntry(h1,"p/f(parabola funciton fit)*f");
    //else leg.AddEntry(h1,"p/f(linear function fit)*f");
    
    // set estimate title
    int stind1 = htitle.Index("_",2);
    int stind2 = htitle.Index("_",stind1+1);
    TString strvar = htitle(stind1+1,stind2-stind1-1);
    if (htitle.Contains("weight2")) leg.SetHeader(strvar+"  Linear");
	else leg.SetHeader(strvar+"  Parabola");
	leg.AddEntry(h1,"Estimation");
    leg.AddEntry(h2,"MC");
    leg.SetBorderSize(0);
    leg.Draw();
    // ratio
    c4->cd(2);
    c4->GetPad(2)->SetGridy();
    c4->GetPad(2)->SetPad(0.0,0.0,1,0.3);
    c4->GetPad(1)->SetPad(0.0,0.3,1,1);
    c4->GetPad(1)->SetTicks();
    c4->GetPad(2)->SetTicks();
    h_copy->GetYaxis()->SetRangeUser(0,2);
    h_copy->GetYaxis()->SetTitle("Estimation/MC");
    h_copy->GetYaxis()->CenterTitle();
    h_copy->SetLineColor(kBlack);
    h_copy->GetXaxis()->SetLabelSize(0);
    h_copy->GetYaxis()->SetTitleSize(0.1);
    h_copy->GetYaxis()->SetLabelSize(0.1);
    h_copy->GetYaxis()->SetTitleOffset(0.5);
    //h_copy->SetMarkerStyle(20);
    // Y = 1 LINE
    static TF1 *f1 = 0;
    if (!f1) f1 = new TF1("f1","1",-1000,1000);
    //f1->SetLineWidth(2);
    h_copy->Draw("E1");
    f1->Draw("same");
    //h_copy->Draw("e1same");
    if (io==1) c4->Print((fileName+"[").data());
    c4->Print(fileName.data());
    if (io==2) c4->Print((fileName+"]").data());
    delete c4;
    delete h_copy;
}