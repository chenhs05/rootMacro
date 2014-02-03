gROOT->Reset();
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>

#include "math.h"
#include "string.h"
#include "tdc_ch_values.cpp"
#include "include/myCali.h"

char *directory = "/home/huangshan/Documents/measurement/SiPM_Measurements/20130611/hv_f3/";
char fileName[1024];

const int chT1 = 8;
const int chTot1 = 12;
const int chT2 = 10;
const int chTot2 = 14;
const float binTDC = 24.414; //ps
//const float binTDC = 25; //ps
const int tWin = 409*10; //10ns*10
const int tTotWin = 818*5;  //20ns*5
const int edgeT[2] = {1,0};
const int edgeTot[2] = {1,0};

const int maxTs = 2097152;//2^21
const int nMax = 1e9;
float t1[nMax],tot1[nMax],t2[nMax],tot2[nMax],tot1t[nMax],tot2t[nMax];
int nEvents;

float tot1Range[2];
float tot2Range[2];
const float tLimit[2]={-2.0,5.0};
const int tBins = (tLimit[1]-tLimit[0])*40;
//const int tBins = (tLimit[1]-tLimit[0])*1000/binTDC;
const float nSigma = 2.0;
const float nTotSigma = 1.0;
const float resTDC = 0.040;

const int NLCorr = 0;

int modeSelect = 1;
int nSelect = 12000;

int noPopUp = 0;

FILE* fscanResultOut;

/************  creat the correction maps for NL correction **********************************
*********************************************************************************************/
void CreatTdcCorrection(char *file,tdc_ch_values *chan,int channel,int edge)
{
	int chIn,tsIn,edgeIn;
	printf("####	Create TDC NL correction for file:\n%s\nchannel: %d	edge: %d\n",file,channel,edge);
	ifstream is;
	is.open(file,ios::in);
	while(!is.eof())
	{
		is>>chIn>>tsIn>>edgeIn;
		if(chIn == channel && edgeIn == edge)	chan->fill(tsIn);
	}
	chan->update_histos();
	is.close();
}

/************  Read .dat file and convert to .root file***************************************
*************  correction with class tdc_ch_value     ****************************************
*********************************************************************************************/
void ReadFile()
{
	printf("\n####	Start to read file with TDC NL correction!	####\n");
	int chIn,tsIn,edgeIn;
	int t1Buf1,t1Buf2,tot1Buf1,tot1Buf2,t2Buf1,t2Buf2,tot2Buf1,tot2Buf2; 
	char buf1[1024];
	char file1[1024];

	sprintf(file1,"%s%s.dat",directory,fileName);

//	tdc_ch_values *chanT1,*chanT2;
//	chanT1 = new tdc_ch_values(1,1024);
//	chanT2 = new tdc_ch_values(1,1024);
//	CreatTdcCorrection(file1,chanT1,chT1,edgeT[0]); 
//	CreatTdcCorrection(file1,chanT2,chT2,edgeT[0]); 
//	chanT1->update_histos();
//	chanT2->update_histos();
//	printf("----> finish create TDC NL correction!	<---- \n");

	cout<<"####	Read File	"<<file1<<"	####"<<"\n";
	ifstream is;
	is.open(file1,ios::in);

	is.clear();
	is.seekg(0,ios::beg);
	nEvents = 0;

	unsigned int i = 0;
	t1Buf1 = 0;
	t1Buf2 = 0;
	tot1Buf1 = 0;
	tot1Buf2 = 0;
	t2Buf1 = 0;
	t2Buf2 = 0;
	tot2Buf1 = 0;
	tot2Buf2 = 0;
	tsBuf = 0;
	int t1Slow = 0;
	int t2Slow = 0;
	
//	FILE* fout;
//	fout=fopen("outData.txt","w");

	while(!is.eof())
	{
		is>>chIn>>tsIn>>edgeIn;
		if(tsIn<tsBuf-1e5)	{
			t1Buf1 = 0;
			t1Buf2 = 0;
			tot1Buf1 = 0;
			tot1Buf2 = 0;
			t2Buf1 = 0;
			t2Buf2 = 0;
			tot2Buf1 = 0;
			tot2Buf2 = 0;
			t1Slow = 0;
			t2Slow = 0;
		}	
		tsBuf=tsIn;
		//if(i%10000 == 0)	cout<<"!!!!    Reading Data	"<<i*10000<<"!!!!	nEvents = "<<nEvents<<"	t1buf2= "<<t1Buf2<<"\n";
		i++;
		switch(chIn)	{
			case chT1:	{
				if(edgeIn == edgeT[0])	{
//					if(NLCorr == 1)	t1Buf1 = chanT1->get_dith_value(tsIn);
//					else	t1Buf1 = tsIn;	
					t1Buf1 = tsIn;

					if(t1Slow == 1 && tot1Buf1 - t1Buf1 > 0 && tot1Buf1-t1Buf1<tTotWin)	{
						t1Slow = 0;
						t1Buf2 = t1Buf1;
					}
					else	{
						t1Slow = 0;
					}
				}
				break;
			}
			case chT2:	{
				if(edgeIn == edgeT[0])	{
//					if(NLCorr == 1)	t2Buf1 = chanT2->get_dith_value(tsIn);
//					else	t2Buf1 = tsIn;
					t2Buf1 = tsIn;
					
					if(t2Slow == 1 && tot2Buf1 - t2Buf1 > 0 && tot2Buf1-t2Buf1<tTotWin)	{
						t2Slow = 0;
						t2Buf2 = t2Buf1;
					}
					else	{
						t2Slow = 0;
					}
				}
				break;
			}
			case chTot1:	{
				if(edgeIn == edgeTot[0])	{
					tot1Buf1 = tsIn;
					if(tot1Buf1>t1Buf1 && abs(tot1Buf1 - t1Buf1) < tTotWin)	{
		      				t1Buf2 = t1Buf1;
					}
					else	{
						t1Slow = 1;
					}
				}
				else if(tot1Buf1 != 0)	{
					tot1Buf2 = tsIn - tot1Buf1;
					if(tot2Buf2 !=0 && abs(tot1Buf1 - tot2Buf1) < tWin && t1Buf2!=0 && t2Buf2!=0)	{
						t1[nEvents] = t1Buf2*binTDC/1000.0;	//ns
						tot1[nEvents] = tot1Buf2*binTDC/1000.0;
						t2[nEvents] = t2Buf2*binTDC/1000.0;
						tot2[nEvents] = tot2Buf2*binTDC/1000.0;
						tot1t[nEvents] = tot1Buf1*binTDC/1000.0;	
						tot2t[nEvents] = tot2Buf1*binTDC/1000.0;	
						nEvents++;

//						fprintf(fout,"%d	%d	%d	%d	%d	%d\n",t1Buf2,t2Buf2,tot1Buf2,tot2Buf2,tot1Buf1,tot2Buf1);
						
						t1Buf1 = 0;
						t1Buf2 = 0;
						tot1Buf1 = 0;
						tot1Buf2 = 0;
						t2Buf1 = 0;
						t2Buf2 = 0;
						tot2Buf1 = 0;
						tot2Buf2 = 0;
						t1Slow = 0;
						t2Slow = 0;
					}
				}
				break;
			}
			case chTot2:	{
				if(edgeIn == edgeTot[0])	{
					tot2Buf1 = tsIn;
					if(tot2Buf1>t2Buf1 && abs(tot2Buf1 - t2Buf1) < tTotWin)	{
		      				t2Buf2 = t2Buf1;
					}
					else	{
						t2Slow = 1;
					}
				}
				else if(tot2Buf1 != 0)	{
					tot2Buf2 = tsIn - tot2Buf1;
					if(tot1Buf2 !=0 && abs(tot1Buf1 - tot2Buf1) < tWin && t1Buf2!=0 && t2Buf2!=0)	{
						t1[nEvents] = t1Buf2*binTDC/1000.0;	//ns
						tot1[nEvents] = tot1Buf2*binTDC/1000.0;
						t2[nEvents] = t2Buf2*binTDC/1000.0;
						tot2[nEvents] = tot2Buf2*binTDC/1000.0;
						tot1t[nEvents] = tot1Buf1*binTDC/1000.0;	
						tot2t[nEvents] = tot2Buf1*binTDC/1000.0;	
						nEvents++;

//						fprintf(fout,"%d	%d	%d	%d	%d	%d\n",t1Buf2,t2Buf2,tot1Buf2,tot2Buf2,tot1Buf1,tot2Buf1);

						t1Buf1 = 0;
						t1Buf2 = 0;
						tot1Buf1 = 0;
						tot1Buf2 = 0;
						t2Buf1 = 0;
						t2Buf2 = 0;
						tot2Buf1 = 0;
						tot2Buf2 = 0;
						t1Slow = 0;
						t2Slow = 0;
					}
				}
				break;
			}
			default:	{
				break;
			}
		}
	}
	is.close();
	
	cout<<"####	Complete reading file :"<<file1<<" !	####"<<"\n";
	cout<<"####	Good Events = "<<nEvents<<", total events = "<<i<<", ratio = "<<1.0*nEvents/i<<"\n";
	cout<<"####	Start to create root file!	####"<<""<<"\n";

	//create root file
	float newT1,newTot1,newT2,newTot2,newTot1T,newTot2T;
	sprintf(buf1,"%s%s.root",directory,fileName);
	TFile *rootFile1 = new TFile(buf1,"RECREATE","Root file from TDC data");
	TTree *h1 = new TTree("h1","Data from TDC");
	h1->Branch("t1",&newT1,"t1/F");
	h1->Branch("tot1",&newTot1,"tot1/F");
	h1->Branch("t2",&newT2,"t2/F");
	h1->Branch("tot2",&newTot2,"tot2/F");
	h1->Branch("tot1t",&newTot1T,"tot1t/F");
	h1->Branch("tot2t",&newTot2T,"tot2t/F");
	for(i=0;i<nEvents;i++)  {
		newT1 = t1[i];
		newTot1 = tot1[i];
		newT2 = t2[i];
		newTot2 = tot2[i];
		newTot1T = tot1t[i];
		newTot2T = tot2t[i];
		h1->Fill();
	}
	//h1->Print();
	h1->Write();
	//fclose (fout);

	//close the root file
	rootFile1->Close();
	cout<<"####	Function ReadFile complete	####"<<"\n\n";
}

/************  Read .root file and plot the spectrum *****************************************
*********************************************************************************************/
void PlotSpectrum()
{
	double mean,sigma,para[16],rang[2],err;
	gStyle->SetOptStat(1111);
	gStyle->SetOptFit(111);
	TFile *f;

	f = new TFile(Form("%s%s.root",directory,fileName));
	TTree *T = (TTree*)gDirectory->Get("h1");
	//T->Print();	

	if(noPopUp)	{
		gROOT->ProcessLine("gROOT->SetBatch()");
	}
	TCanvas *c1 = new TCanvas("CTR measurement","CTR measurement results",40,20,1000,700);
	c1->Divide(3,2);
	c1->cd(1);
	T->Draw(Form("tot1>>TOT1(%d,%f,%f)",1000,0.0,1000.0));
	c1->cd(2);
	T->Draw(Form("tot2>>TOT2(%d,%f,%f)",1000,0.0,1000.0));
	TSpectrum *s = new TSpectrum(2);
	int nfound = s->Search(TOT1,2,"",0.05);
	float *xpeak = s->GetPositionX();
	TOT1->Fit("gaus","NQ");
	gaus->SetRange(xpeak[0]-50,xpeak[0]+50);
	TOT1->Fit("gaus","RQ");
	gaus->GetParameters(&para[0]);
	tot1Range[0]=para[1]-nTotSigma*para[2];
	tot1Range[1]=para[1]+nTotSigma*para[2];
	TCut totcut1 = Form("tot1>%.1f && tot1<%.1f ",tot1Range[0],tot1Range[1]);
	//TCut totcut1 = Form("tot1>%.1f && tot1<%.1f ",para[1]-nTotSigma*para[2],para[1]+nTotSigma*para[2]);

	nfound = s->Search(TOT2,2,"",0.05);
	xpeak = s->GetPositionX();
	TOT2->Fit("gaus","NQ");
	gaus->SetRange(xpeak[0]-50,xpeak[0]+50);
	TOT2->Fit("gaus","RQ");
	gaus->GetParameters(&para[0]);
	tot2Range[0]=para[1]-nTotSigma*para[2];
	tot2Range[1]=para[1]+nTotSigma*para[2];
	TCut totcut2 = Form("tot2>%.1f && tot2<%.1f ",tot2Range[0],tot2Range[1]);
	//TCut totcut2 = Form("tot2>%.1f && tot2<%.1f ",para[1]-nTotSigma*para[2],para[1]+nTotSigma*para[2]);

	c1->cd(3);
	T->Draw(Form("tot1>>TOT1_aftre_cut"),totcut1);
	TOT1_aftre_cut->Rebin(5);
	c1->cd(4);
	T->Draw(Form("tot2>>TOT2_aftre_cut"),totcut2);
	TOT2_aftre_cut->Rebin(5);
	c1->cd(5);
//	TF1 *f1 = new TF1("f1","gaus");
	if(modeSelect ==1)	{
		TCut totcut = totcut1 && totcut2 && Form("Entry$ < %d", nSelect);
	}
	else	{
		TCut totcut = totcut1 && totcut2;	
	}
	T->Draw(Form("t1-t2>>CTR1(%d,%f,%f)",tBins,tLimit[0],tLimit[1]),totcut);
	CTR1->Fit("gaus","NQ");
//	gaus->GetParameters(&para[0]);
//	gaus->SetRange(para[1]-nSigma*para[2], para[1]+nSigma*para[2]);
  	mean=CTR1->GetMean();
  	sigma=CTR1->GetRMS();
	gaus->SetRange(mean-nSigma*sigma,mean+nSigma*sigma);
	gaus->GetRange(rang[0],rang[1]);
	cout<<"range:	"<<rang[0]<<"	"<<rang[1]<<endl;
	CTR1->Fit("gaus","NQR");
	gaus->GetParameters(&para[0]);
	gaus->SetRange(para[1]-nSigma*para[2], para[1]+nSigma*para[2]);
	gaus->GetRange(rang[0],rang[1]);
	cout<<"range:	"<<rang[0]<<"	"<<rang[1]<<endl;
	CTR1->Fit("gaus","NQR");
	gaus->GetParameters(&para[0]);
	gaus->SetRange(para[1]-nSigma*para[2], para[1]+nSigma*para[2]);
	gaus->GetRange(rang[0],rang[1]);
	cout<<"range:	"<<rang[0]<<"	"<<rang[1]<<endl;
	CTR1->Fit("gaus","R");
	CTR1->GetXaxis()->SetTitle("delay (ns)");
	CTR1->GetYaxis()->SetTitle("Counts");
	CTR1->GetXaxis()->SetLabelSize(0.05);
	CTR1->GetYaxis()->SetLabelSize(0.05);
	CTR1->SetTitle("CTR");

	gaus->GetParameters(&para[0]);
	err=gaus->GetParError(2);
	
	TF1 *ftemp = new TF1("ftemp","[0]*exp(-0.5*((x-[1])/[2])**2)",tLimit[0],tLimit[1]);
	ftemp->SetParameters(para[0],para[1],para[2]);
	ftemp->SetLineWidth(3);
	ftemp->SetLineStyle(2);
	ftemp->SetLineColor(2);
	ftemp->Draw("same");
	c1->Update();
	///////////save histogram to .root file
	TFile *fHist = new TFile(Form("%s%s-hist.root",directory,fileName),"recreate");
	c1->Write();
	fHist->ls();
	fHist->Close();
	//////////
	
	double ctrRes,ctrErr;
	ctrRes=sqrt(para[2]*para[2]-resTDC*resTDC)*2.35*1000.0;
	ctrErr=err*2.35*1000.0;
	fprintf(fscanResultOut,"%f	%f\n",ctrRes,ctrErr);

	//////////////output fit result to txt
//	ofstream outResult;
//	outResult.open(Form("%s%s-fitResult.txt",directory,fileName),ios::app);
//	outResult<<"########    "<<fileName<<"    #######"<<"\n";
//	outResult<<"nSigma	NLCorr	tLimit[0]	tLimit[1]	tBins	tot1Range[0]	tot1Range[1]	tot2Range[0]	tot2Range[1]"<<"\n";
//	outResult<<nSigma<<'	'<<NLCorr<<'	'<<tLimit[0]<<'	'<<tLimit[1]<<'	'<<tBins<<'	'<<tot1Range[0]<<'	'<<tot1Range[1]<<'	'<<tot2Range[0]<<'	'<<tot2Range[1]<<"\n";
//	outResult<<"Sigma[ps]	Err_of_Sigma[ps]"<<"\n";
//  	outResult<<para[2]*1000.0<<"\t\t"<<err*1000.0<<"\n\n\n";
//	outResult.close();
//	/////////////////////////
	
	//wait to close the file
	//wait();
	
	f->Close();	
	c1->Close();
	if(noPopUp)	{
		gROOT->ProcessLine("gROOT->SetBatch(kFALSE)");
	}
}

/************ wait funtion !!  *************************************************************** 
*********************************************************************************************/
void wait(){

	TTimer *timer = new TTimer("gSystem->ProcessEvents();", 50, kFALSE);
	char *input;
	Bool_t done = kFALSE;
	do{
		timer->TurnOn();
		timer->Reset();
		input=Getline("Type <return> to continue : ");
		timer->TurnOff();
		if(input){			done =kTRUE;		}	
	}while(!done);
}

void DoCali(){
	int i,j;
	int nCount,nEvent;
	const int nTotMax=1e8;
	float ct[nTotMax],tot[2][nTotMax],newT1,newTot1,newT2,newTot2,newTot1T,newTot2T;
	char psName[1024];
	TFile *f;

	f = new TFile(Form("%s%s.root",directory,fileName));
	TTree *T = (TTree*)gDirectory->Get("h1");
	nEvent= h1->GetEntries();
	h1->SetBranchAddress("t1",&newT1);
	h1->SetBranchAddress("tot1",&newTot1);
	h1->SetBranchAddress("t2",&newT2);
	h1->SetBranchAddress("tot2",&newTot2);
	h1->SetBranchAddress("tot1t",&newTot1T);
	h1->SetBranchAddress("tot2t",&newTot2T);

	for(i=0;i<nEvent+2;i++){
		ct[i] = 0;
		tot[0][i] = 0;
		tot[1][i] = 0;
	}

	nCount = 0;
	if(modeSelect ==1 ) nEvent = nSelect;
	for(i=1; i<nEvent; i++){
		h1->GetEntry(i);
		if(newTot1>tot1Range[0] && newTot1<tot1Range[1] && newTot2>tot2Range[0] && newTot2<tot2Range[1]){
			ct[nCount] = newT1-newT2;
			tot[0][nCount] = newTot1;
			tot[1][nCount] = newTot2;
			nCount++;
		}
	}

	sprintf(psName,"CTR");
	int nCharge = 2;
	int nCor = 5;
	float tLimit = 5.0;
	int tBinsCali = 2*tLimit*40;
	float minTot[2],maxTot[2],fitMin[2],fitMax[2];
	minTot[0] = tot1Range[0];
	minTot[1] = tot2Range[0];
	maxTot[0] = tot1Range[1];
	maxTot[1] = tot2Range[1];
	fitMin[0] = tot1Range[0];
	fitMin[1] = tot2Range[0];
	fitMax[0] = tot1Range[1];
	fitMax[1] = tot2Range[1];
	int totBins = 20;
	//TF1 *fitf = new TF1("fitf","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6+[7]*x**7+[8]*x**8+[9]*x**9+[10]*x**10");
	TF1 *fitf = new TF1("fitf","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6");
	DoCorrection(ct,tot[0],nCount,nCharge,nCor,psName,tLimit,tBinsCali,minTot,maxTot,totBins,fitMin,fitMax,fitf,nSigma,nTotMax);
}

/************ Main funtion !!  *************************************************************** 
*********************************************************************************************/
void tdc_ctr_scan()
{
	cout<<"\n\n";
	cout<<"#############################################################################\n";
	cout<<"####################	Start of tdc_ctr analysis	#####################\n";
	cout<<"#############################################################################\n";
	
	char buf[1024];
	sprintf(buf,"%sscanResult.txt",directory);
	fscanResultOut = fopen(buf,"w");

	remove(Form("%scaliResult.txt",directory));

	float hv;
	int i,j;
	for(i=1;i<=16;i+=1)	{
		//for(j=59;j<=83;j+=2)	{
			hv=i*0.2+71.3;
			//sprintf(fileName,"chip1_bias_%d",i);
			sprintf(fileName,"chip2_%d_%.1f",i,hv);
			sprintf(buf,"%s%s.dat",directory,fileName);
			ifstream isTest;
			isTest.open(buf,ios::in);
			if(!isTest) {
				cout<<"\n"<<buf<<" doesn't exist!!"<<"\n"<<"\n";
				//return;
				continue;
			}
			else {
				isTest.close();
				cout<<"File:"<<"\n"<<"########  "<<buf<<"  ########"<<"\n\n";
				fprintf(fscanResultOut,"%.1f	",hv);
				//fprintf(fscanResultOut,"%d	",i);
				ReadFile();
				PlotSpectrum();
				DoCali();
			}
		//}
	}
	fclose(fscanResultOut);
	//free(fileName);

	cout<<"#############################################################################\n";
	cout<<"####################	End of tdc_ctr analysis		####################\n";
	cout<<"#############################################################################\n";
	cout<<"\n\n";
}
