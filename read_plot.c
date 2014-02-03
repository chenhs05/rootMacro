gROOT->Reset();
#include <Riostream.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include "math.h"
#include "string.h"
#include "TROOT.h"
#include "Riostream.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"  
#include "TRandom.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPostScript.h"

const Int_t nMax=1000000;
Float_t var[nMax],w[nMax],var1,varbuf,;

Char_t *directory="/home/huangshan/Documents/_measurement/SiPM_Measurements/20121022WithSTiC";
Char_t *fileName="La.Ba.Th.35.5.12";
Char_t buf1[1024],buf2[1024];

Float_t varUpt,varLow,varMax,varMin;
Int_t nVar;
const Int_t nBins=1000;


void read_plot()
{
	Int_t i,j,k;
//	for(i=0;i<nMax;i++)	w[i]=1;
	
	sprintf(buf1,"%s\/%s.txt",directory,fileName);
	ifstream indata;
	indata.open(buf1,ios::in);
	
	if(!indata) {
    cout<<endl<<buf1<<" doesn't exist!!"<<endl<<endl;
		exit(-1);
  }
  else {
  	sprintf(buf2,"%s\/%s.root",directory,fileName);	//root
		TFile *file1 = new TFile(buf2,"RECREATE");
		TTree *h1 = new TTree("h1","data from File");
		h1->Branch("var",&var1,"var/F");
		
		varMax=0;
		varMin=0;
		nVar=0;
		
		while(!indata.eof())  {
			indata>>varbuf>>var1;
//			indata>>var1;
			
//			if(var1>0.132) continue;
			
			h1->Fill();
			
			//var[nVar]=var1*1e3;
			var[nVar]=var1*1e12;
			if(varMax<var[nVar])  varMax=var[nVar];
			if(varMin>var[nVar])  varMin=var[nVar];
				
			nVar++;
			if(nVar==nMax)	break;
		}
		h1->Write();
		indata.close();
		file1->Close();
		
		varbuf=varMax-varMin;
		varMax=varMax+0.025*varbuf;
		varMin=varMin-0.025*varbuf;
		
		TH1D *spec = new TH1D("SiPM","Spectrum",nBins,varMin,varMax);
		//TH1D *spec = new TH1D("SiPM","Spectrum",nBins*0.05,71000,76000);
		for(i=0;i<nVar;i++)  {
			spec->Fill(var[i]);
		}
		
		gStyle->SetOptFit(1111);
		gStyle->SetOptStat(1111);
		TCanvas *c1 = new TCanvas("c1","c1",10,10,950,750);
		spec->Draw();
		//spec->Fit("gaus");
		spec->GetXaxis()->SetTitle("t [ps]");
		//spec->GetXaxis()->SetTitle("Maximum of SiPM signals [mV]");
		spec->GetXaxis()->CenterTitle();
		spec->GetYaxis()->SetTitle("Counts");
		spec->GetYaxis()->CenterTitle();
  }
	
}
