#!/bin/python
import csv
from ROOT import TGraph, TCanvas, TFile, TString
from ROOT import TLegend, TMultiGraph, TH1F, TGraphAsymmErrors
from array import array
def getEffi(fileName):
    arr1 = array('f',[])
    arr2 = array('f',[])
    with open(fileName) as fff:
        row = csv.reader(fff,delimiter='\t')
        for num in row:
            arr1.append(float(num[0]))
            arr2.append(float(num[1]))

    return [arr1,arr2]


def main():
    c1 = TCanvas("c1","c1",3)
    histList = array('f',[500,700,900,1100,1300,1500,1800,2200,2700])
    passh = TH1F('passh','passh',8,histList);
    allh = TH1F('allh','allh',8,histList);
    tgas = TGraphAsymmErrors()
    zpList = array('f',[600,800,1000,1200,1400,1700,2000,2500])
    effiList = array('f',[])
    pList = array('f',[])
    allList = array('f',[])
    FourJetCorr = []
    FiveJetCorr = []
    method = ['min #DeltaR (A0#rightarrowbb)','min #DeltaMh','min pt As (A0#rightarrowbb)','max pt As (zp#rightarrowh+A0)','max #DeltaR (zp#rightarrowh+A0)']
    for i in range(5):
        FourJetCorr.append(array('f',[]))
        FiveJetCorr.append(array('f',[]))

    for mzp in zpList:
        fName = "corr_MZp" + str(int(mzp)) + "_MA0300.txt"
        arr=getEffi(fName)
        # correctness
        for i in range(5):
            FourJetCorr[i].append(arr[0][i+1]/arr[0][0])
            FiveJetCorr[i].append(arr[1][i+1]/arr[1][0])
        
        # effi
        effiList.append(arr[0][6]/10000)
        pList.append(arr[0][6])
        allList.append(10000)
        ## effi error bar
        passh.Fill(mzp,arr[0][6])
        allh.Fill(mzp,10000)
    
    # effi plot
    tgas.BayesDivide(passh,allh)
    tgas.GetXaxis().SetTitle('MZp (GeV)')
    tgas.GetYaxis().SetTitle('Efficiency')
    tgas.Draw('ALP')
    c1.Print('effi.pdf')
    ## tgEffi = TGraph(8,zpList,effiList)
    ## tgEffi.Draw('ALP')
    
    # corr plot
    tgList4 = [TGraph(8,zpList,FourJetCorr[i]) for i in range(5) ]
    tgList5 = [TGraph(8,zpList,FiveJetCorr[i]) for i in range(5) ]
    mg = TMultiGraph()
    mg.SetMaximum(1.)
    lg = TLegend(0.59,0.55,0.9,0.9)
    lg.AddEntry(0,'4 Leading Jets','')
    for i in range(5):
        tgList5[i].SetLineStyle(9)
        tgList4[i].SetLineWidth(2)
        tgList5[i].SetLineWidth(2)
        tgList4[i].SetLineColor(51+8*i)
        tgList5[i].SetLineColor(51+8*i)
        mg.Add(tgList4[i])
        mg.Add(tgList5[i])
        lg.AddEntry(tgList4[i],method[i],'lp')
    lg.AddEntry(0,'5 Leading Jets','')
    for i in range(5):
        lg.AddEntry(tgList5[i],method[i],'lp')

    c1.Print('Correctness.pdf[')
    mg.Draw('ALP')
    mg.GetXaxis().SetTitle('MZp (GeV)')
    mg.GetYaxis().SetTitle('Correctness')
    lg.Draw()
    c1.Print('Correctness.pdf')
    c1.Print('Correctness.pdf]')

if __name__=="__main__":
    main()
