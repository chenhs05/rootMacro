gROOT->Reset();
#include <Riostream.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>

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

//Char_t *directory = "/home/huangshan/Documents/measurement/SiPM_Measurements/data_ana";
//Char_t *fileName="20121115-1";

char *directory = "/home/huangshan/Documents/measurement/SiPM_Measurements/20121203TDC_Scope/Scope";
char *fileName = "20121203";

const Int_t nMax=500000;
Float_t t[nMax],tot[nMax],data_buf[64];
Int_t nEvents;
Char_t buf[1024],buf2[1024],psName[1024];

const Float_t minT = 69.5;
const Float_t maxT = 75.5;
const Float_t nBinsAllT = 200;
const Float_t minTot = 20;
const Float_t maxTot = 160;
const Float_t nBinsAllTot = 140;

const Float_t shift_t = -35.6;
const Int_t nCor = 5;
const Float_t nSigma = 1.5; 
const Float_t halfPeakWidth = 1;
/*//20121105-1
const Int_t nBinsPeakT = 200;
const Float_t minPeakTot[2] = {43,0};
const Float_t maxPeakTot[2] = {57,0};
const Int_t nBinsPeakTot = 20;
const Float_t fitMinTot[2] = {43,0};
const Float_t fitMaxTot[2] = {57,0};
*/
////20121105-2
//const Int_t nBinsPeakT = 200;
//const Float_t minPeakTot[2] = {47,0};
//const Float_t maxPeakTot[2] = {63,0};
//const Int_t nBinsPeakTot = 20;
//const Float_t fitMinTot[2] = {47,0};
//const Float_t fitMaxTot[2] = {63,0};
//20121115-1
const Int_t nBinsPeakT = 59;
const Float_t minPeakTot[2] = {52,0};
const Float_t maxPeakTot[2] = {65,0};
const Int_t nBinsPeakTot = 20;
const Float_t fitMinTot[2] = {minPeakTot[0],0};
const Float_t fitMaxTot[2] = {maxPeakTot[0],0};

Int_t nCharge = 1;


//TF1 *fitf = new TF1("fitf","[0]+[1]/sqrt(x)+[2]/x+[3]*x");  //��������
//TF1 *fitf = new TF1("fitf","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6");
TF1 *fitf = new TF1("fitf","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6+[7]*x**7+[8]*x**8+[9]*x**9+[10]*x**10");
/******************************************************************************************************************************
//�������� by Huangshan Chen, ��д/��װ�� �������� wjb_test.c
//ttm��ʱ�����ݵ��׵�ַ��һά����Ļ�����������
//tot��������ݵ��׵�ַ����ά���ݵĻ���tot[0];
//nCount: ��Чevent�ĸ����� 
//nCharge����Ҫ�����ĵ�ɸ�����
//nCor������������
//psName�����ps�ļ����ļ����������������directory Ŀ¼������psName_Calibration.ps�ļ���
//tLimit��float��ʽ��ʱ��ֲ�ͼ�ķ�Χ�� [-tLimit,tLimit]��
//tBins: ʱ��ֲ���bin����
//totMin: ���Ȼ���time over threshold�ֲ�����Сֵ�����ڶ�T-Qͼ��x��ķ�Χ����Ӧ��ͬ�ĵ�ɵķ��ȷ�Χ����������׵�ַ��
//totMin: ���Ȼ���time over threshold�ֲ������ֵ�����ڶ�T-Qͼ��x��ķ�Χ����������׵�ַ��
//totBins: T-Q��ϵͼx���bin����
//fitMin�� �����ķ�Χ����Сֵ����������׵�ַ��
//fitMax�� �����ķ�Χ�����ֵ����������׵�ַ��
//fitf�� ����������
//Nsigma�� ��ʱ��ֲ�����Gaussion fit�ķ�Χ��ȡƽ��ֵ����Nsigma��sigma
//nTotMax��ttm����tot[0]�ĳ��ȣ�����ȡ����ȷ��totֵ����Ϊtot[1][0]��tot[0][0]�ĵ�ַ���nTotMax
//IsSplit���Ƿ�RPC�ķ��ȷֳ�avalanche �� streamer �ֱ�����ٷ���һ��fit�õ�ʱ��ֱ��ʣ�0--��1--�ǣ�Ĭ��Ϊ0���������Ҫ������
//totSplit��avalanche��streamer���ȵķֽ�㣬Ĭ��Ϊ0���������Ҫ��������
//IdxRpcCharge��RPC�ĵ����tot��ά�����ϵ���ţ���RPC��ɴ�tot[IdxRpcCharge][0]��ʼ��Ĭ��Ϊ4���������Ҫ��������
******************************************************************************************************************************/
		 DoCorrection(t2,tot2,nPeakCounts,nCharge,nCor,psName,halfPeakWidth,nBinsPeakT,minPeakTot,maxPeakTot,nBinsPeakTot,fitMinTot,fitMaxTot,fitf,nSigma,nPeakCounts);
void DoCorrection(Float_t *ttm,Float_t *tot,Int_t nCount,Int_t nCharge,Int_t nCor,Char_t *psName,Float_t tLimit,Int_t tBins,Float_t *totMin,Float_t *totMax,Int_t totBins,Float_t *fitMin,Float_t *fitMax,TF1 *fitf,Float_t Nsigma,Int_t nTotMax,Int_t IsSplit=0,Float_t totSplit=0,Int_t IdxRpcCharge=4)
{
	cout<<"######## "<<psName<<" Calibration ########"<<endl;
	cout<<"Good events:    "<<nCount<<endl;
	
	Double_t mean,sigma,para[16],para2[16],err;
  TH1D *T;
  TH2D *Corr;
  TH1D *htemp;
  Int_t nPad;
  Double_t cor;
  Char_t buf[1024];
  
	TF1 *f1 = new TF1("f1","gaus");
	TF1 *ftemp = new TF1("ftemp","[0]*exp(-0.5*((x-[1])/[2])**2)",-tLimit[0],tLimit[0]);
	f1->SetLineColor(2);
	ftemp->SetLineWidth(3);
	ftemp->SetLineStyle(2);
	ftemp->SetLineColor(2);
	
	sprintf(buf,"%s\/%s_Calibration.ps",directory,psName);
	TPostScript *ps = new TPostScript(buf,112);
  ps->Range(24,16);
  
  TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,800);
	c1->SetFillColor(kWhite);
  c1->GetFrame()->SetFillColor(21);
 	c1->GetFrame()->SetBorderSize(6);
  c1->GetFrame()->SetBorderMode(-1);
  
	//gStyle->SetOptFit(110);
	//gStyle->SetOptStat(110);
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	
	//////////////output fit result
  ofstream outResult;
  sprintf(buf,"%s\/%s_FitResult.txt",directory,fileName);
  outResult.open(buf,ios::app);
  outResult<<"########    "<<fileName<<"    #######"<<endl;
  outResult<<"nSigma	nCor	nBinsPeakT	nBinsPeakTot	minPeakTot[0]	maxPeakTot[0]	fitMinTot[0]	fitMaxTot[0]"<<endl;
  outResult<<nSigma<<'	'<<nCor<<'	'<<nBinsPeakT<<'	'<<nBinsPeakTot<<'	'<<minPeakTot[0]<<'	'<<maxPeakTot[0]<<'	'<<fitMinTot[0]<<'	'<<fitMaxTot[0]<<endl;
  outResult<<"nCorr	Sigma	Err of Sigma"<<endl;
  /////////////////////////
  
  for(Int_t i=0;i<nCor;i++)
  {
  	cout<<"########### The No."<<i+1<<" correction. ##########"<<endl;
  	ps->NewPage();
  	c1->Clear();
  	c1->Divide(1,1);
  	c1->cd(1)->Update();
  	//c1->cd(1)->SetLogy();
  	
  	sprintf(buf,"%s_%d_Original",psName,i);
//		sprintf(buf," ");
  	T= new TH1D(buf,buf,tBins,-tLimit,tLimit); 
	  T->SetLineWidth(2);
	  
  	for(Int_t m=0;m<nCount;m++)
  	{
  		T->Fill(ttm[m]);
  	}
  	mean=T->GetMean();
  	sigma=T->GetRMS();
	  f1->SetRange(mean-Nsigma*sigma, mean+Nsigma*sigma);
	  T->Fit("f1","NQR");
	  f1->GetParameters(&para[0]);
	  f1->SetRange(para[1]-Nsigma*para[2], para[1]+Nsigma*para[2]);
	  T->Fit("f1","NQR");
	  f1->GetParameters(&para[0]);
	  f1->SetRange(para[1]-Nsigma*para[2], para[1]+Nsigma*para[2]);
	  T->Draw();
	  T->Fit("f1","QR");
	  f1->GetParameters(&para[0]);
		err=f1->GetParError(2);
	  
	  ftemp->SetParameters(para[0],para[1],para[2]);
	  ftemp->Draw("same");
	  
//	  T->GetXaxis()->SetTitle("((PMT1+PMT2)-(PMT3+PMT4))/4 [��35ps]");
		T->GetXaxis()->SetTitle("t [ns]");
		T->GetXaxis()->CenterTitle();
		T->GetYaxis()->SetTitle("Counts");
		T->GetYaxis()->CenterTitle();
		/////////////////////
		if(i==0)
			{
  			outResult<<"No.0"<<'\t'<<para[2]*1000.0<<'\t'<<err*1000.0<<endl;
  		}
  	////////////////////
		
		c1->cd(1)->Update();
		
  	for(Int_t j=0;j<nCharge;j++)
  	{
  		if(j%2==0) 
  		{
  			ps->NewPage();
  			c1->Clear();
  			c1->Divide(2,2);
  			nPad=1;
  		}
  		else 
  		{
  			nPad=2;
  		}
  		sprintf(buf,"%s_Correction(%d)_%d",psName,i,j);
  		Corr = new TH2D(buf,buf,totBins,totMin[j],totMax[j],tBins,-tLimit,tLimit);
  		
  		for(m=0;m<nCount;m++)
  		{
  			Corr->Fill(*(tot+j*nTotMax+m),ttm[m]);
  		}
  		c1->cd(nPad);
  		Corr->Draw("colz");
  			
  		Corr->GetXaxis()->SetTitle("TOT [ns]");
			Corr->GetXaxis()->CenterTitle();
			Corr->GetYaxis()->SetTitle("t [ns]");
			Corr->GetYaxis()->CenterTitle();
			c1->cd(nPad)->Update();
  		
  		Corr->ProfileX();
			sprintf(buf,"%s_Correction(%d)_%d_pfx",psName,i,j);
			htemp=(TProfile *)gDirectory->Get(buf);
			c1->cd(nPad)->Update();
//			htemp->Draw("QRsame");
//			htemp->SetLineColor(kBlack);
//			htemp->SetMarkerColor(kBlack);
//			htemp->SetMarkerStyle(21);
//			htemp->SetMarkerSize(0.6);
//  		c1->cd(nPad)->Update();
  		
  		if(IsSplit==1&&j==IdxRpcCharge)
  		{
  			TF1 *fitf2= new TF1("fitf2","fitf");
	  		htemp->Fit("fitf","q","QRsame",fitMin[j],totSplit);
	  		fitf->GetParameters(&para[0]);
	  		fitf2->SetRange(fitMin[j],totSplit);
	  		fitf2->SetParameters(&para[0]);
				fitf2->SetLineWidth(2);
	  		fitf2->Draw("same");
  			c1->cd(nPad)->Update();
	  		htemp->Fit("fitf","q","QRsame",totSplit,fitMax[j]);
	  		fitf->GetParameters(&para2[0]);
  			c1->cd(nPad)->Update();	  		
	  		for(m=0;m<nCount;m++)
	  		{
	  			if((*(tot+j*nTotMax+m))<totSplit)
	  			{
//		  			cor=para[0]+para[1]/sqrt((*(tot+j*nTotMax+m)))+para[2]/(*(tot+j*nTotMax+m))+para[3]*(*(tot+j*nTotMax+m));
						cor=0;
		  			ttm[m]=ttm[m]-cor;
	  			}
	  			if((*(tot+j*nTotMax+m))>totSplit)
	  			{
//		  			cor=para2[0]+para2[1]/sqrt((*(tot+j*nTotMax+m)))+para2[2]/(*(tot+j*nTotMax+m))+para2[3]*(*(tot+j*nTotMax+m));
						cor=0;
		  			ttm[m]=ttm[m]-cor;
	  			}
	  		} 				
  		}
  		else
  		{
	  		//htemp->Fit("fitf","q","QRsame",fitMin[j],fitMax[j]);
	  		htemp->Fit("fitf","Nq","QRsame",fitMin[j],fitMax[j]);
	  		c1->cd(nPad)->Update();
	  		fitf->GetParameters(&para[0]);
	  		for(m=0;m<nCount;m++)
	  		{
//	  			cor=para[0]+para[1]/sqrt((*(tot+j*nTotMax+m)))+para[2]/(*(tot+j*nTotMax+m))+para[3]*(*(tot+j*nTotMax+m));
	  			cor=fitf->Eval((*(tot+j*nTotMax+m)));
	  			ttm[m]=ttm[m]-cor;
	  		}
  		}
  		
//  		sprintf(buf,"%s-T(%d)_%d",psName,i,j);
  		sprintf(buf," ");
  		T= new TH1D(buf,buf,tBins,-tLimit,tLimit); 
  		for(m=0;m<nCount;m++)
  		{
  			T->Fill(ttm[m]);
  		}
  		c1->cd(nPad+2);
//  		c1->cd(nPad+2)->SetLogy();
			c1->cd(nPad+2)->SetBorderSize(6);
	  	mean=T->GetMean();
	  	sigma=T->GetRMS();
		  f1->SetRange(mean-Nsigma*sigma, mean+Nsigma*sigma);
  		T->Fit("f1","NQR");
		  f1->GetParameters(&para[0]);
		  f1->SetRange(para[1]-Nsigma*para[2], para[1]+Nsigma*para[2]);
		  T->Fit("f1","NQR");
		  f1->GetParameters(&para[0]);
		  f1->SetRange(para[1]-Nsigma*para[2], para[1]+Nsigma*para[2]);
		  T->Draw();
		  T->Fit("f1","QR");
		  f1->GetParameters(&para[0]);
		  err=f1->GetParError(2);
		  
		  ftemp->SetParameters(para[0],para[1],para[2]);
		  ftemp->Draw("same");
		  
		  T->GetXaxis()->SetTitle("t [ns]");
//		  T->GetXaxis()->SetTitle("t [\\times 35ps]");
//		  T->GetXaxis()->SetTitle("((PMT1+PMT2)-(PMT3+PMT4))/4 [\\times 35ps]");
			T->GetXaxis()->CenterTitle();
			T->GetYaxis()->SetTitle("Counts");
			T->GetYaxis()->CenterTitle();
  		c1->cd(nPad+2)->Update();
  		
  		/////////////////////
  		outResult<<"No."<<i+1<<'\t'<<para[2]*1000.0<<'\t'<<err*1000.0<<endl;
  		////////////////////
  	}
  }
  
  outResult.close();
	ps->Close();
	delete c1;
}

void ReadFile(Char_t *file1)
{
	cout<<"####	Read File	"<<buf<<"	####"<<endl;
	ifstream is;
	is.open(file1,ios::in);
	nEvents=0;
	while(is>>data_buf[0]>>data_buf[1])
	{
		t[nEvents]=data_buf[0];
		tot[nEvents]=data_buf[1];
		nEvents++;
//		cout<<nEvents<<endl;
	}
	is.close();
	cout<<"####	Complete reading file "<<buf<<" !	####"<<endl<<endl;
}

void PlotStatus()
{
	cout<<"####	Plot the spectrums	####"<<endl;
	int i,j,k;
	TObjArray histlist(0);
	TH1D *histT = new TH1D("t","t",nBinsAllT,minT,maxT);
	TH1D *histTOT = new TH1D("TOT","TOT",nBinsAllTot,minTot,maxTot);
	histlist.Add(histT);
	histlist.Add(histTOT);
	
//	cout<<endl<<"nEvents = "<<nEvents<<endl<<endl;
	for(i=0;i<nEvents;i++)
	{
		histT->Fill(t[i]);
		histTOT->Fill(tot[i]);
//		cout<<t[i]<<'\t'<<tot[i]<<endl;
	}
	
	TCanvas *cStatus1 = new TCanvas("status1","status1",10,10,800,600);
//  cStatus->Clear();
//	cStatus->Divide(2,1);
//	cStatus->cd(1);
	histT->Draw();
	histT->GetXaxis()->SetTitle("delay [ns]");
	histT->GetXaxis()->CenterTitle();
	histT->GetYaxis()->SetTitle("Counts");
	histT->GetYaxis()->CenterTitle();
	cStatus1->Update();
	
	TCanvas *cStatus2 = new TCanvas("status2","status2",10,10,800,600);
//	cStatus->cd(2);
	histTOT->Draw();
	histTOT->GetXaxis()->SetTitle("TOT [ns]");
	histTOT->GetXaxis()->CenterTitle();
	histTOT->GetYaxis()->SetTitle("Counts");
	histTOT->GetYaxis()->CenterTitle();
	cStatus2->Update();
	
	sprintf(buf2,"%s\\%s-RunStatus.root",directory,fileName);
	TFile *hf = new TFile(buf2,"recreate");
	histlist.Write();
	cout<<"####	Complete ploting spectrums !	####"<<endl<<endl;
}

void TQCorrection()
{
	cout<<"####	Start T-Q correction!	####"<<endl;
	Int_t i;
	Int_t nPeakCounts=0;
	Float_t t2[nMax],tot2[nMax];
	for(i=0;i<nEvents;i++)
	{
//		cout<<tot[i]<<endl;
		if(tot[i]>minPeakTot[0] && tot[i]<maxPeakTot[0])
			{
				t2[nPeakCounts]=t[i]+shift_t;
				tot2[nPeakCounts]=tot[i];
				nPeakCounts++;
			}
	}
	cout<<"	Number of good events = "<<nPeakCounts<<endl;
	sprintf(psName,"Scope");
	DoCorrection(t2,tot2,nPeakCounts,nCharge,nCor,psName,halfPeakWidth,nBinsPeakT,minPeakTot,maxPeakTot,nBinsPeakTot,fitMinTot,fitMaxTot,fitf,nSigma,nPeakCounts);
	cout<<"####	Complete T-Q correction!	####"<<endl<<endl;
}

void cali_osci()
{
	sprintf(buf,"%s\/%s.txt",directory,fileName);	//�ļ���·��
	/*
	TFile *file1 = new TFile(buf);
	if(file1->IsZombie()) {
    cout<<endl<<buf<<" doesn't exist!!"<<endl<<endl;
//    exit(-1);
  }*/
  ifstream isTest;
	isTest.open(buf,ios::in);
	if(!isTest) {
    cout<<endl<<buf<<" doesn't exist!!"<<endl<<endl;
    return;
  }
  else {
  	isTest.close();
  	cout<<"File:"<<endl<<"########  "<<buf<<"  ########"<<endl<<endl;
  	ReadFile(buf);
		PlotStatus();
		TQCorrection();
	}
}
