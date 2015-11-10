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
	fileName.push_back("/home/huangshan/gitRepos/ipython_notebook/201510_t_vs_tot_result_0.txt");

	ifstream infile;
	int nEntr = fileName.size();
	for(int i=0;i<nEntr;++i)
	{
		TGraphErrors *g;

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
			int style_count=0;
			double max_pulse_amp=1.000;
			n = 0;
			g = new TGraphErrors();
			infile.seekg(0,ios::beg);
			while(true)	{
				infile>>var[0]>>var[1]>>var[2];
				infile>>var[3]>>var[4]>>var[5];
				infile>>var[6]>>var[7];
				infile>>var[8]>>var[9];
				infile>>var[10];

				if(infile.eof()) break;				// end of file, break

				if(var[10] > 5) {
					x=var[1]*330; // injected charge in pC
					y=var[2];
					xErr=0;
					yErr=var[3];
				
					g->SetPoint(n,x,y);
					g->SetPointError(n,xErr,yErr);
					n++;
					cout << n << endl;
				}

				if(var[1] == max_pulse_amp) {
					g->SetMarkerSize(1.0);
					g->SetMarkerStyle(20+style_count);
					g->SetMarkerColor(2+style_count);
					g->SetLineColor(2+style_count);

					graph_list.push_back(g);
					g = new TGraphErrors();
					style_count++;
					n=0;
				}

			}
		}
		infile.close();

	}


	TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
	TMultiGraph *mg = new TMultiGraph();
//	c1->SetGrid();

	nEntr=graph_list.size();
	cout << "Number of plots: " << nEntr << endl;
	for(i=0;i<nEntr;++i)
	{
		mg->Add(graph_list.back());
		graph_list.pop_back();
	}
	mg->Print();
	mg->Draw("ap");
	mg->SetTitle("ToT vs. C");
	mg->GetXaxis()->SetTitle("C_{in}[pC]");
	mg->GetYaxis()->SetTitle("ToT [s]");
//	mg->GetXaxis()->CenterTitle();
//	mg->GetYaxis()->CenterTitle();
	c1->BuildLegend();
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
