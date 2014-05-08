//root macro for sptr analysis with HPTDC
//create by H.Chen
//#include <TROOT.h>
//gROOT->Reset();
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

#include "TTree.h"
#include "TGraph.h"
#include "TFile.h"
#include "TLine.h"

char *directory = "/home/huangshan/Documents/measurement/SiPM_Measurements/201401_sptr/goodEvent/";
char *fileName = "sptr_17_2_1";

const int chLaser = 1;
const int chStic = 4;
const int edgeLaser[2] = {0,1};
const int edgeStic[2] = {1,0};

const float binTDC = 25.0; ///24.414; //ps
const int totHyst = 409; //
const float maxDelay = 4096;	//*25ps

const int nMax = 1e7; 
float t[nMax],tot[nMax],totLaser[nMax];
int nEvents;

const float minT = 35;
const float maxT = 45;
const float binningAllT = 0.05;
int nBinsAllT = int((maxT-minT)/binningAllT);
const float minTot = 0;
const float maxTot = 200;
const float binningAllTot = 0.5;
int nBinsAllTot = int((maxTot-minTot)/binningAllTot);

const float shift_t = 40;
const int nCor = 5;
const float nSigma = 1.5; 
const float halfPeakWidth = 10;

float binningPeakT = 0.05;
int nBinsPeakT = int(2*halfPeakWidth/binningPeakT);
double minPeakTot[2] = {56,0};
double maxPeakTot[2] = {64,0};
float binningPeakTot = 0.5;
int nBinsPeakTot = int((maxPeakTot[0]-minPeakTot[0])/binningPeakTot);
//const int nBinsPeakTot = 20;
double fitMinTot[2] = {minPeakTot[0]+1,0};
double fitMaxTot[2] = {maxPeakTot[0]-1,0};

int nCharge = 1;

int nSelect = 13000;
bool modeSelect = 0;

bool NLCorr = 0;
/************  Read .dat file and convert to .root file, without NL correction *************
*********************************************************************************************/
/*void ReadFile(char *file1)
{
	cout<<"############## nMax = "<<nMax<<"\n";
	int chIn,tsIn,edgeIn,dt;
	float tsBuf,tsLaserBuf,totLaserBuf,tsSticBuf,tSticBuf,tsSticBuf2,totSticBuf;
	char buf1[1024];
	float tStic,totStic,totTrigger;
	cout<<"####	Read File	"<<file1<<"	####"<<"\n";
	ifstream is;
	is.open(file1,ios::in);
	nEvents = 0;

	tsBuf=0;
	tsLaserBuf = 0;
	tsSticBuf  = 0;
	tSticBuf   = 0;
	tsSticBuf2 = 0;
	totSticBuf = 0;

	int i = 0;
	while(is>>chIn>>tsIn>>edgeIn)
	{
		if(tsIn<tsBuf)	{
			tsLaserBuf = 0;
			tsSticBuf  = 0;
			tSticBuf   = 0;
			tsSticBuf2 = 0;
			totSticBuf = 0;
		}	
		tsBuf=tsIn;
		if(i%100000 == 0)	cout<<"!!!!    Reading Data	"<<i<<"!!!!	nEvents = "<<nEvents<<"\n";
		i++;
		switch(chIn)	{
			case chLaser:	{
				if(edgeIn == edgeLaser[0])	tsLaserBuf = tsIn;
				break;
			}
			case chStic:		{
				if(edgeIn == edgeStic[0])	{
					if(totSticBuf == 0)	{
						tsSticBuf = tsIn;
						tSticBuf = tsSticBuf - tsLaserBuf;
					}
					else	{	
						dt = tsIn-tsSticBuf2;
						if(dt>totHyst)	{
							if(tSticBuf<maxDelay)	{
								tot[nEvents] = totSticBuf*binTDC/1000.0;	//ns
								t[nEvents] = tSticBuf*binTDC/1000.0;	//ns
								nEvents++;
							}
							totSticBuf = 0;
							tsSticBuf = tsIn;
							tSticBuf = tsSticBuf - tsLaserBuf;
						}
					}
				}
				else if(tsSticBuf != 0)	{
					totSticBuf = tsIn-tsSticBuf;
					tsSticBuf2 = tsIn;
				}	
				break;		
			}
			default:	{
		//		cout<<chIn<<"	!!!!!	switch channel == default	!!!!!!!!!"<<"\n";
				break;
			}
		}	
	}
	is.close();
	
	cout<<"####	Complete reading file :"<<file1<<" !	####"<<"\n";
	cout<<"####	Start to create root file!	####"<<""<<"\n";

	//create root file
	sprintf(buf1,"%s\/%s.root",directory,fileName);
	TFile *rootFile1 = new TFile(buf1,"RECREATE","Root file from TDC data");
	TTree *h1 = new TTree("h1","Data from TDC");
	h1->Branch("t",&tStic,"t/F");
	h1->Branch("tot",&totStic,"tot/F");
	for(i=0;i<nEvents;i++)  {
		tStic = t[i];
		totStic = tot[i];
		h1->Fill();
	}
	h1->Print();
	h1->Write();

	//close the root file
	rootFile1->Close();
	cout<<"####	Function ReadFile complete	####"<<"\n";
}*/

/*********** create the NL correction mapping for HPTDC ************************************
 ********************************************************************************************/
void CreatTdcCorrection(char *file,tdc_ch_values *chan,int channel,int edge)
{
	int chIn,tsIn,edgeIn;
	printf("####	Create TDC NL correction for file:\n%s\n",file);
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
	long int chIn,tsIn,edgeIn,dt;
	long int tsBuf,tsLaserBuf,totLaserBuf,tsSticBuf,tSticBuf,tsSticBuf2,totSticBuf,tsLaserBuf2;
	char buf1[1024];
	char file1[1024];
	float tStic,totStic,totTrigger;
	sprintf(file1,"%s%s.dat",directory,fileName);

	tdc_ch_values *chanLaser, *chanStic;	
	if(NLCorr == 1)	{
		printf("\n####	Start to read file with TDC NL correction!	####\n");
		chanLaser = new tdc_ch_values(1,1024);
		chanStic = new tdc_ch_values(1,1024);
		CreatTdcCorrection(file1,chanLaser,chLaser,edgeLaser[0]); 
		CreatTdcCorrection(file1,chanStic,chStic,edgeStic[0]); 
		chanLaser->update_histos();
		chanStic->update_histos();
		printf("----> finish create TDC NL correction!	<---- \n");
	}

	cout<<"####	Read File	"<<file1<<"	####"<<"\n";
	ifstream is;
	is.open(file1,ios::in);

	is.clear();
	is.seekg(0,ios::beg);
	nEvents = 0;

	tsBuf=0;
	tsLaserBuf = 0;
	tsLaserBuf2 = 0;
	tsSticBuf  = 0;
	tSticBuf   = 0;
	tsSticBuf2 = 0;
	totSticBuf = 0;
	int i = 0;
	while(!is.eof())
	{
		is>>chIn>>tsIn>>edgeIn;
		if(tsIn<tsBuf-1e5)	{
			tsLaserBuf = 0;
			tsLaserBuf2 = 0;
			tsSticBuf  = 0;
			tSticBuf   = 0;
			tsSticBuf2 = 0;
			totSticBuf = 0;
		}	
		tsBuf=tsIn;
		if(i%100000 == 0)	cout<<"!!!!    Reading Data	"<<i<<"!!!!	nEvents = "<<nEvents<<"\n";
		i++;
		switch(chIn)	{
			case chLaser:	{
				if(edgeIn == edgeLaser[0]){
					if(NLCorr == 1) tsLaserBuf = chanLaser->get_dith_value(tsIn);
					else tsLaserBuf = tsIn;
				//	tsSticBuf  = 0;
				//	tsSticBuf2 = 0;
				//	tSticBuf = 0;
				//	totSticBuf = 0;

				}
				break;
			}
			case chStic:		{
				if(edgeIn == edgeStic[0])	{
					if(totSticBuf == 0)	{
					if(tsLaserBuf != tsLaserBuf2)	{
						if(NLCorr == 1) tsSticBuf = chanStic->get_dith_value(tsIn);
						else tsSticBuf = tsIn;
						if(tsSticBuf > tsLaserBuf )	{
							tSticBuf = tsSticBuf - tsLaserBuf;
							tsLaserBuf2 = tsLaserBuf;
						}
						else tsSticBuf = 0;
					}
					}
					else	{	
						if(NLCorr ==1 ) dt = chanStic->get_dith_value(tsIn)-tsSticBuf2;
						else dt = tsIn-tsSticBuf2;
						if(dt>totHyst)	{
							if(tSticBuf<maxDelay)	{
								tot[nEvents] = totSticBuf*binTDC*0.001;	//ns
								t[nEvents] = tSticBuf*binTDC*0.001;	//ns

					//			if(t[nEvents]< 30)	printf("%f %f %d %d %d\n",t[nEvents],tot[nEvents],tsSticBuf,tsLaserBuf2,tsLaserBuf);

								nEvents++;
							}
							totSticBuf = 0;
							tsSticBuf = 0;
							if(tsLaserBuf != tsLaserBuf2)	
							{
								if(NLCorr == 1) tsSticBuf = chanStic->get_dith_value(tsIn);
								else tsSticBuf = tsIn;
								if(tsSticBuf > tsLaserBuf )	{
									tSticBuf = tsSticBuf - tsLaserBuf;
									tsLaserBuf2 = tsLaserBuf;
								}
								else tsSticBuf = 0;

							}
						}
					}
				}
				else if(tsSticBuf != 0)	{
					if(NLCorr == 1) {
						totSticBuf = chanStic->get_dith_value(tsIn)-tsSticBuf;
						tsSticBuf2 = chanStic->get_dith_value(tsIn);
					} 
					else {
						totSticBuf = tsIn-tsSticBuf;
						tsSticBuf2 = tsIn;
					}
				}	
				break;		
			}
			default:	{
		//		cout<<chIn<<"	!!!!!	switch channel == default	!!!!!!!!!"<<"\n";
				break;
			}
		}	
	}
	is.close();
	
	cout<<"####	Complete reading file :"<<file1<<" !	####"<<"\n";
	cout<<"####	Start to create root file!	####"<<""<<"\n";

	//create root file
	sprintf(buf1,"%s%s.root",directory,fileName);
	TFile *rootFile1 = new TFile(buf1,"RECREATE","Root file from TDC data");
	TTree *h1 = new TTree("h1","Data from TDC");
	h1->Branch("t",&tStic,"t/F");
	h1->Branch("tot",&totStic,"tot/F");
	for(i=0;i<nEvents;i++)  {
		tStic = t[i];
		totStic = tot[i];
		h1->Fill();
	}
	h1->Print();
	h1->Write();

	//close the root file
	rootFile1->Close();
	cout<<"####	Function ReadFile complete	####"<<"\n";
}

/************  Read .root file and put the data into variable*********************************
*********************************************************************************************/
void ReadRootFile(TFile *rootFile,char *treeName)
{
	float tBuf,totBuf;
 	TTree *h1 = (TTree*)rootFile->Get(treeName);
	h1->Print();
   	h1->SetBranchAddress("t",&tBuf);
   	h1->SetBranchAddress("tot",&totBuf);
	if(modeSelect == 1)	nEvents = nSelect;
	else	nEvents = (int)h1->GetEntries();
	for(int i=0;i<nEvents;i++)
	{
		h1->GetEntry(i);
		t[i]	= tBuf;
		tot[i]	= totBuf;
	}
	rootFile->Close();
}

/************ Plot the status of the data!  ************************************************** 
*********************************************************************************************/
void PlotStatus()
{
	cout<<"\n"<<"####	Plot the spectrums	####"<<"\n";
	int i,j,k;
	char buf2[1024];
//	TObjArray histlist(0);
	TH1D *histT = new TH1D("t","t",nBinsAllT,minT,maxT);
	TH1D *histTOT = new TH1D("TOT","TOT",nBinsAllTot,minTot,maxTot);
	TH2D *tVsTot = new TH2D("tVsTot","t vs. TOT",nBinsPeakTot,minPeakTot[0],maxPeakTot[0],nBinsPeakT,shift_t-halfPeakWidth,shift_t+halfPeakWidth);
	TGraph *sigmaSlice; 
//	histlist.Add(histT);
//	histlist.Add(histTOT);
//	histlist.Add(sigmaSlice);
//	histlist.Add(tVsTot);
	
//	cout<<"\n"<<"nEvents = "<<nEvents<<"\n"<<"\n";
	for(i=0;i<nEvents;i++)
	{
		histT->Fill(t[i]);
		histTOT->Fill(tot[i]);
//		cout<<t[i]<<'\t'<<tot[i]<<"\n";
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
	double maxY;
	maxY = gPad->GetUymax(); 
//	cout<<"max of Y axis = "<<maxY<<"\n";
	TLine *l1 = new TLine(minPeakTot[0],0,minPeakTot[0],maxY);
	l1->SetLineColor(kRed);
	l1->Draw();
	TLine *l2 = new TLine(maxPeakTot[0],0,maxPeakTot[0],maxY);
	l2->SetLineColor(kRed);
	l2->Draw();
	cStatus2->Update();
	

//	TF1 *f1 = new TF1("f1","gaus");
//	const int maxNSlice = 1024;
//	TH1D *histSlice[maxNSlice];
//	double para[16],mean,sigma,err;
//	double binSliceLow[maxNSlice],binSliceUp[maxNSlice],meanBinSlice[maxNSlice],meanT[maxNSlice],sigmaT[maxNSlice],sigmaSigmaT[maxNSlice],nEventsSlice[maxNSlice];
//	for(i=0;i<nBinsPeakTot;i++)	{
//		sprintf(buf2,"tSlice_%d",i);
//		histSlice[i] = new TH1D(buf2,buf2,int((maxT-minT)/binningPeakT),minT,maxT);
//		binSliceLow[i]=minPeakTot[0]+i*binningPeakTot;
//		binSliceUp[i]=minPeakTot[0]+(i+1)*binningPeakTot;
//		meanT[i]=0;
//		sigmaT[i]=0;
//		sigmaSigmaT[i]=0;
//		nEventsSlice[i]=0;
//	}
//	for(i=0;i<nEvents;i++)	{
////		if(i%100000==0)	cout<<"!!!!!!!!!!!!!	"<<i<<"	!!!!!!!!!!!!!!!!!!!"<<"\n";
//		for(int j=0;j<nBinsPeakTot;j++)	{
//			if(tot[i]>=binSliceLow[j] && tot[i]<binSliceUp[j])	{
//				histSlice[j]->Fill(t[i]);
//				nEventsSlice[j]++;
//			}
//		}
//	}
//	ofstream outProfile;
//	sprintf(buf2,"%s%s_SliceProfile.txt",directory,fileName);
//	outProfile.open(buf2,ios::trunc);
//	outProfile<<"Slice#	TOTLow	TOTUP	nEvents	meanT	sigmaT	sigmaSigmaT"<<"\n";
//	for(i=0;i<nBinsPeakTot;i++)	{
//		mean=histSlice[i]->GetMean();
//		sigma=histSlice[i]->GetRMS();
//		f1->SetRange(mean-nSigma*sigma, mean+nSigma*sigma);
//		histSlice[i]->Fit("f1","NQR");
//		f1->GetParameters(&para[0]);
//		f1->SetRange(para[1]-nSigma*para[2], para[1]+nSigma*para[2]);
//		histSlice[i]->Fit("f1","NQR");
//		f1->GetParameters(&para[0]);
//		f1->SetRange(para[1]-nSigma*para[2], para[1]+nSigma*para[2]);
//		histSlice[i]->Fit("f1","NQR");
//		f1->GetParameters(&para[0]);
//		err=f1->GetParError(2);
//		meanT[i]=para[1];
//		sigmaT[i]=para[2]*1000.0;
//		sigmaSigmaT[i]=err*1000.0;
//		meanBinSlice[i]=0.5*(binSliceLow[i]+binSliceUp[i]);
//		outProfile<<i<<"	"<<binSliceLow[i]<<"	"<<binSliceUp[i]<<"	"<<nEventsSlice[i]<<"	"<<para[1]<<"	"<<para[2]*1000.0<<"	"<<err*1000.0<<"\n";
//	}
//	outProfile.close();
//	TCanvas *cStatus3 = new TCanvas("status3","status3",10,10,800,600);
//	sigmaSlice = new TGraph(nBinsPeakTot,meanBinSlice,sigmaT);
//	sigmaSlice->SetTitle("sigmaT vs tot");
//	sigmaSlice->Draw("A*");
//	cStatus3->Update();

	for(i=0;i<nEvents;i++)	{
		if(tot[i]>minPeakTot[0] && tot[i]<maxPeakTot[0])	{	
			tVsTot->Fill(tot[i],t[i]);
		}
	}
	TCanvas *cStatus4 = new TCanvas("status4","status4",10,10,800,600);
	tVsTot->Draw("colz");
  	tVsTot->GetXaxis()->SetTitle("TOT [ns]");
	tVsTot->GetXaxis()->CenterTitle();
	tVsTot->GetYaxis()->SetTitle("t [ns]");
	tVsTot->GetYaxis()->CenterTitle();

//	TProfile *tTotProfile = tVsTot->ProfileX();
//	tTotProfile->Draw("QRsame");
//	tTotProfile->SetMarkerStyle(21);
//	tTotProfile->SetMarkerSize(0.6);
	cStatus4->Update();

	sprintf(buf2,"%s%s_RunStatus.root",directory,fileName);
	TFile *hf = new TFile(buf2,"recreate");
//	histlist.Write();
	cStatus1->Write();
	cStatus2->Write();
//	cStatus3->Write();
	cStatus4->Write();

	hf->Close();
	cout<<"####	Complete ploting spectrums !	####"<<"\n"<<"\n";
}

/************  Selecte the data and do the  calibration. *************************************
*********************************************************************************************/
void TQCorrection()
{
	cout<<"####	Start T-Q correction!	####"<<"\n";
	int i;
	int nPeakCounts=0;
	double t2[nMax],tot2[nMax];
	char psName[1024];
	sprintf(psName,"%s",fileName);
	
	for(i=0;i<nEvents;i++)
	{
//		cout<<tot[i]<<"\n";
		if(tot[i]>minPeakTot[0] && tot[i]<maxPeakTot[0])
			{
				t2[nPeakCounts]=t[i]-shift_t;
				tot2[nPeakCounts]=tot[i];
				nPeakCounts++;
			}
	}

	cout<<"	Number of good events = "<<nPeakCounts<<"\n";
	sprintf(psName,"TDC");
	TF1 *fitf = new TF1("fitf","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6+[7]*x**7+[8]*x**8+[9]*x**9+[10]*x**10");
	DoCorrection(t2,tot2,nPeakCounts,nCharge,nCor,psName,halfPeakWidth,nBinsPeakT,minPeakTot,maxPeakTot,nBinsPeakTot,fitMinTot,fitMaxTot,fitf,nSigma,nPeakCounts);
	cout<<"####	Complete T-Q correction!	####"<<"\n"<<"\n";
}

/************ Main funtion !!  *************************************************************** 
*********************************************************************************************/
void tdc_sptr()
{
	char buf[1024];
	char modeIn;
	cout<<"*******************************************************************"<<"\n";
	cout<<"*******************************************************************"<<"\n";
	cout<<"The mode of this macro:"<<"\n";
	cout<<"1. Read .dat file and convert to .root file."<<"\n";
	cout<<"2. Read .root file and plot the t and TOT spectrum."<<"\n";
	cout<<"3. Read .root file and do TOT-t calibration."<<"\n";
	cout<<"4. Read .dat file, convert to .root file and do TOT-t calibration."<<"\n";
	cout<<"*******************************************************************"<<"\n";
	cout<<"Please select the mode: ";
	cin>>modeIn;
	cout<<"\n";
	switch(modeIn)	{
		case '1':	{
			sprintf(buf,"%s%s.dat",directory,fileName);	
		
			ifstream isTest;
			isTest.open(buf,ios::in);
			if(!isTest) {
				cout<<"\n"<<buf<<" doesn't exist!!"<<"\n"<<"\n";
				return;
			}
			else {
				isTest.close();
				cout<<"File:"<<"\n"<<"########  "<<buf<<"  ########"<<"\n"<<"\n";
				ReadFile();
				PlotStatus();
			}
			break;
		}
		case '2':	{
			sprintf(buf,"%s%s.root",directory,fileName);	
		
			TFile *f1 = new TFile(buf);
			if(f1->IsZombie())	{
				cout<<"Error openning file!!"<<"\n";
				return;
			}
			else {
				cout<<"File:"<<"\n"<<"########  "<<buf<<"  ########"<<"\n"<<"\n";
				ReadRootFile(f1,"h1");
				PlotStatus();
			}
			break;
		}
		case '3':	{
			sprintf(buf,"%s%s.root",directory,fileName);	
		
			TFile *f1 = new TFile(buf);
			if(f1->IsZombie())	{
				cout<<"Error openning file!!"<<"\n";
				return;
			}
			else {
				cout<<"File:"<<"\n"<<"########  "<<buf<<"  ########"<<"\n"<<"\n";
				ReadRootFile(f1,"h1");
			//	PlotStatus();
				TQCorrection();
				f1->Close();
			}
			break;
		}
		case '4':	{
			sprintf(buf,"%s%s.dat",directory,fileName);	
			ifstream isTest;
			isTest.open(buf,ios::in);
			if(!isTest) {
				cout<<"\n"<<buf<<" doesn't exist!!"<<"\n"<<"\n";
				return;
			}
			else {
				isTest.close();
				cout<<"File:"<<"\n"<<"########  "<<buf<<"  ########"<<"\n"<<"\n";
				ReadFile();
				PlotStatus();
				TQCorrection();
			}
			break;
		}
		default:	{
			cout<<"Wrong Mode!!!!"<<"\n";
			break;
		}
	}
}
