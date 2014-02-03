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
const Int_t nVar=2;
Float_t var[nVar][nMax];

Char_t *directory="";
Char_t *fileName1="F6_00001";
Char_t *fileName2="F5_00001";
Char_t buf1[1024],buf2[1024];

Float_t varUp[nVar],varLow[nVar],varMax[nVar],varMin[nVar],var1,varbuf;
Int_t nEntrs,indx1,indxMin[nVar],indxMax[nVar],indxStart[nVar];
const Int_t nBins=1000;


void read_plot()
{
	Int_t i,j,k;
//	for(i=0;i<nMax;i++)	w[i]=1;
	
	ifstream indata[nVar];

	sprintf(buf1,"%s\\%s.txt",directory,fileName1);
	indata[0].open(buf1,ios::in);
	sprintf(buf2,"%s\\%s.txt",directory,fileName2);
	indata[1].open(buf2,ios::in);
	
	if(!indata[0]) {
    	cout<<endl<<buf1<<" doesn't exist!!"<<endl<<endl;
		exit(-1);
  	}
	else if(!indat[1])  {
	cout<<endl<<buf2<<" doesn't exist!!"<<endl<<endl;
		exit(-1);
  	}
  	else {
  		/*sprintf(buf2,"%s\\%s.root",directory,fileName);	//root
		TFile *file1 = new TFile(buf2,"RECREATE");
		TTree *h1 = new TTree("h1","data from File");
		h1->Branch("var",&var1,"var/F");*/
		
		

		Int_t nCount;		
		for(i=0;i<nVar;i++)  {
			varMax[i]=0;
			varMin[i]=nMax;

			indxMin[i]=nMax;
			indxMax[i]=0;
			nCount=0;
			while(!indata[0].eof())  {
				indata>>indx1>>var1;
	//			indata>>var1;

	//			if(var1>0.132) continue;
				//h1->Fill();
	//			var[nVar]=var1*1e3;

				if(indxMin[i]>indx1)	indxMin[i]=indx1;
				if(indxMax[i]<indx1)	indxMax[i]=indx1;
				var[i][nCount]=var1*1e12;
				if(varMax[i]<var[i][nCount])  varMax[i]=var[i][nCount];
				if(varMin[i]>var[i][nCount])  varMin[i]=var[i][nCount];
				
				nCount++;
				if(nCount==nMax)	break;
			}
			//h1->Write();
			indata[i].close();
			//file1->Close();
		
			varbuf=varMax[i]-varMin[i];
			varMax[i]=varMax[i]+0.025*varbuf;
			varMin[i]=varMin[i]-0.025*varbuf;
		}

		indxStart[0]=indxMin[0]<=indxMin[1]?0:indxMin[0]-indxMin[1];
		indxStart[1]=indxMin[1]<=indxMin[0]?0:indxMin[1]-indxMin[0];
		
		TH1D *spec = new TH1D("SiPM","Spectrum",nBins,varMin,varMax);
//		TH1D *spec = new TH1D("SiPM","Spectrum",nBins*0.05,81000,84000);
		for(i=0;i<nVar;i++)  {
			spec->Fill(var[i]);
		}
		
		gStyle->SetOptFit(1111);
		gStyle->SetOptStat(1111);
		TCanvas *c1 = new TCanvas("c1","c1",10,10,950,750);
		spec->Draw();
//		spec->Fit("gaus");
		spec->GetXaxis()->SetTitle("t [ps]");
//		spec->GetXaxis()->SetTitle("Maximum of SiPM signals [mV]");
		spec->GetXaxis()->CenterTitle();
		spec->GetYaxis()->SetTitle("Counts");
		spec->GetYaxis()->CenterTitle();
  }
	
}
