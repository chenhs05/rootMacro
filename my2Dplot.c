gROOT->Reset();
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>
#include"myWait.c"

char *fileName = "/home/huangshan/Documents/measurement/SiPM_Measurements/20130609/biasScan_2/scanResult_1.0sigma_1.5Esigma.txt";

void my2Dplot()	{
	int n;
	double ctr[1024],err[1024],bias1[1024],bias2[1024];
	ifstream infile;
	infile.open(fileName,ios::in);

	if(!infile) {
		cout<<"\n"<<fileName<<" doesn't exist!!"<<"\n"<<"\n";
		return;
	}
	else {
		n=0;
		while(!infile.eof())	{
			infile>>bias1[n]>>bias2[n]>>ctr[n]>>err[n];
			n++;
		}
	}
	infile.close();

	TGraph2D *g = new TGraph2D("g","CTR;bias1;bias2;CTR",n-1, bias1, bias2, ctr);
	TCanvas *c1 = new TCanvas("c1","c1",0,0,600,400);
	g->Draw("colz");

	wait();
	c1->Close();
}
