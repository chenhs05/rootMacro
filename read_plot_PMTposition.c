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

//PMT 和 RPC的时间，幅度
Float_t Pm1_t,Pm2_t,Pm3_t,Pm4_t,Pm5_t,Pm6_t;
Float_t q3_1,q3_2,q3_3,q3_4,q1,q2,q3,q4;
Float_t Trpc_t[12],Trpc_q[12];

Char_t *directory="E:\\root\\swap";
Char_t *fileName="run619";
Char_t buf[1024];

Float_t min_qpair[2]={800,150}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4")
Float_t max_qpair[2]={2100,600};//
Float_t min_qdiff[2]={-1200,-450};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4")
Float_t max_qdiff[2]={1000,100};//
Float_t min_diff[2]={-95,-20}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t")
Float_t max_diff[2]={70,110}; //
Float_t min_position=-180; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")
Float_t max_position=20;

TObjArray histlist(0);
TH1D *sct1 = new TH1D("(Pm1_t-Pm2_t)/2","Scintillator",100,-50,50);
TH1D *sct2 = new TH1D("(Pm3_t-Pm4_t)/2","Scintillator",100,-50,50);

//用于判断事例是否在所设的PMT cut之内
bool InsideCut()
{
	bool isCut=0;
	if(Pm1_t>0 && Pm2_t>0 && Pm3_t>0 &&Pm4_t>0 
		&& q3_1>=10 && q3_2>=10 && q3_3>=10 && q3_4>=10
		&& q3_1+q3_2 >= min_qpair[0] && q3_1+q3_2 <= max_qpair[0]
		&& q3_3+q3_4 >= min_qpair[1] && q3_3+q3_4 <= max_qpair[1]
		&& q3_1-q3_2 >= min_qdiff[0] && q3_1-q3_2 <= max_qdiff[0]
		&& q3_3-q3_4 >= min_qdiff[1] && q3_3-q3_4 <= max_qdiff[1]
		&& Pm1_t-Pm2_t >= min_diff[0] && Pm1_t-Pm2_t <= max_diff[0]
		&& Pm3_t-Pm4_t >= min_diff[1] && Pm3_t-Pm4_t <= max_diff[1]
		&& (Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)>= min_position && (Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)<= max_position
//		&& (Pm1_t+Pm2_t-Pm3_t-Pm4_t)>min_tof && (Pm1_t+Pm2_t-Pm3_t-Pm4_t)<max_tof
		)		{isCut=1;}
		return isCut;
}

void Read(TFile *f,Char_t *treeName)
{
 	TTree *h1 = f->Get(treeName);
	
	for(Int_t i=0;i<12;i++)
	{
		sprintf(buf,"Trpc_t%d",i+1);
		h1->SetBranchAddress(buf,&Trpc_t[i]);
		sprintf(buf,"Trpc_q%d",i+1);
		h1->SetBranchAddress(buf,&Trpc_q[i]);
	}
   	
   	h1->SetBranchAddress("q3_1",&q3_1);
   	h1->SetBranchAddress("q3_2",&q3_2);
   	h1->SetBranchAddress("q3_3",&q3_3);
   	h1->SetBranchAddress("q3_4",&q3_4);
   	
   	h1->SetBranchAddress("Pm1_t",&Pm1_t);
   	h1->SetBranchAddress("Pm2_t",&Pm2_t);
   	h1->SetBranchAddress("Pm3_t",&Pm3_t);
   	h1->SetBranchAddress("Pm4_t",&Pm4_t);
   	h1->SetBranchAddress("Pm5_t",&Pm5_t);
   	h1->SetBranchAddress("Pm6_t",&Pm6_t);
}

void plotHis()
{
	Int_t i,j,k;
	Int_t nCount=0;
	Int_t nEvent= (Int_t)h1->GetEntries();
	
	for(i=0;i<nEvent;i++)
	{
		h1->GetEntry(i);
		if(InsideCut()==0) continue;
		if(Pm5_t>0&&Pm6_t>0)
			{
				sct1->Fill((Pm1_t-Pm2_t)/2);
				sct2->Fill((Pm3_t-Pm4_t)/2);
			}
	}	
	
	gStyle->SetOptFit(1111);
	gStyle->SetOptStat(1111);
	TCanvas *c1 = new TCanvas("c1","dt1",10,10,1150,750);
	sct1->Draw();
	sct1->Fit("gaus");
	
	TCanvas *c2 = new TCanvas("c2","dt2",10,10,1150,750);
	sct2->Draw();
	sct2->Fit("gaus");
}

void read_plot_PMTposition()
{
	histlist.Add(sct1);
	histlist.Add(sct2);
	
	sprintf(buf,"%s\\%s.root",directory,fileName);	//root文件的路径
	TFile *file1 = new TFile(buf);
	
	if(file1->IsZombie()) {
    cout<<endl<<buf<<" doesn't exist!!"<<endl<<endl;
    exit(-1);
  }
  else {
  	cout<<"Root File:"<<endl<<"########  "<<buf<<"  ########"<<endl<<endl;
  	
  	sprintf(buf,"h1");  //tree的名称
		Read(file1,buf);		
		Int_t nEvents= (Int_t)h1->GetEntries();
		
		plotHis();
  }
}