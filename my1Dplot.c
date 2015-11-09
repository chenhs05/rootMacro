gROOT->Reset();
#include <fstream>
#include <iostream>
#include <iomanip>
#include <list>
//#include <vector>

//#include <stdio.h>
//#include <stdlib.h>
//#include <unistd.h>

void my1Dplot()	{
	std::list<char *> fileName;
	std::list<TGraphErrors *> graph_list;
	int n,buf;
	double var[16];

//	fileName.push_back("/home/huangshan/Dropbox/measurement/stic3_measurements/201407_scan_single_crystal_temp_18deg/20140704_tthr_scan_02_04/result/cspect_fit_32bins.log");
//	fileName.push_back("/home/huangshan/Dropbox/measurement/stic3_measurements/201407_scan_single_crystal_temp_18deg/20140704_tthr_scan_02_04/result/cspect_fit_64bins.log");
//	fileName.push_back("/home/huangshan/Dropbox/measurement/stic3_measurements/201407_scan_single_crystal_temp_18deg/20140704_tthr_scan_02_04/result/cspect_fit_128bins.log");
	fileName.push_back("/home/huangshan/Dropbox/measurement/stic3_measurements/20140517_thr_scan/scan_result/com_spi_3_t_thr_scan.txt");

	ifstream infile;
	int nEntr = fileName.size();
	for(int i=0;i<nEntr;++i)
	{
		TGraphErrors *g = new TGraphErrors();

		//for the list 
		char *last_file_name = fileName.back();
		fileName.pop_back();

		infile.open(last_file_name,ios::in);
		if(!infile) {
			cout<<"\n"<<last_file_name<<" doesn't exist!!"<<"\n"<<"\n";
			return;
		}
		else {
			float x,y,xErr,yErr;
			n=0;
			infile.seekg(0,ios::beg);
			while(!infile.eof())	{
				infile>>var[0]>>var[1]>>var[2];
				infile>>var[3]>>var[4]>>var[5];
				infile>>var[6]>>var[7];
				infile>>var[8]>>var[9];
				infile>>var[10]>>var[11];

				x=181.418-3.19416*var[1]; // threshold dac value to threshold level
				y=var[5];
				xErr=0;
				yErr=var[6];
				
				g->SetPoint(n,x,y);
				g->SetPointError(n,xErr,yErr);
				n++;
			}
		}
		infile.close();

		g->SetMarkerSize(1.0);
		g->SetMarkerStyle(20+i);
		g->SetMarkerColor(2+i);
		g->SetLineColor(2+i);

		graph_list.push_back(g);
	}


	TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
	TMultiGraph *mg = new TMultiGraph();
	c1->SetGrid();

	nEntr=graph_list.size();
	for(i=0;i<nEntr;++i)
	{
		mg->Add(graph_list.back());
		graph_list.pop_back();
	}
	mg->Draw("ap");
	mg->SetTitle("TThr_scan");
	mg->GetXaxis()->SetTitle("TThreshold [fC]");
	mg->GetYaxis()->SetTitle("CTR FWHM [ps]");
//	mg->GetXaxis()->CenterTitle();
//	mg->GetYaxis()->CenterTitle();
	c1->Update();

	wait();
	c1->Close();
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
