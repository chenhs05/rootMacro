gROOT->Reset();
#include <fstream>
#include <iostream>
#include <iomanip>
#include <list>
//#include <vector>

//#include <stdio.h>
//#include <stdlib.h>
//#include <unistd.h>

bool DEBUG = 0;

void my1Dplot()	{
	const EColor colors[] = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink};
	std::list<char *> * fileName = new std::list<char *>;
	std::list<TGraphErrors *> * graph_list = new std::list<TGraphErrors *>;
	int n,buf;
	double var[16];

//	fileName->push_back("/home/huangshan/Dropbox/measurement/stic3_measurements/201407_scan_single_crystal_temp_18deg/20140704_tthr_scan_02_04/result/cspect_fit_32bins.log");
	//fileName->push_back("/home/huangshan/gitRepos/ipython_notebook/result/20151110_t_vs_tot_conv_remove_low_ethr.txt");
	//fileName->push_back("/home/huangshan/gitRepos/ipython_notebook/result/20151110_t_vs_tot_no_low_pass_conv.txt");
	//fileName->push_back("/home/huangshan/gitRepos/ipython_notebook/result/20151110_t_vs_tot_no_low_pass_filter_low_ethr_conv.txt");
	//fileName->push_back("/home/huangshan/gitRepos/ipython_notebook/result/20151115_t_vs_tot_C33p_w_low_pass_filter_conv_remove_entries.txt");
	fileName->push_back("/home/huangshan/gitRepos/ipython_notebook/result/20151111_t_vs_tot_C33p_no_low_pass_filter_conv.txt");

	ifstream infile;
	int nEntr = fileName->size();

	TLegend* leg = new TLegend(0.1, 0.6, 0.3, 0.9);
	char * low_pass;

	for(int i=0;i<nEntr;i++)
	{
		TGraphErrors *g;

		//for the list
		char *last_file_name = fileName->back();
		fileName->pop_back();

		infile.open(last_file_name,ios::in);
		if(!infile) {
			cout<<"\n"<<last_file_name<<" doesn't exist!!"<<"\n"<<"\n";
			return;
		}
		else {
			// read in the first line to get the initial value of last_dac and last_ethr
			infile.seekg(0,ios::beg);
			infile>>var[0]>>var[1]>>var[2];
			infile>>var[3]>>var[4]>>var[5];
			infile>>var[6]>>var[7];
			infile>>var[8]>>var[9];
			infile>>var[10];
			infile>>var[11];
			if(infile.eof()) break;				// end of file, break

			float x,y,xErr,yErr;
			int style_count=0;
			float capa_in=33;
			float last_dac=var[0];
			float last_ethr=var[11];

			int style_diff_count = 0;


			n = 0;
			g = new TGraphErrors();
			infile.seekg(0,ios::beg);
			while(true)	{
				infile>>var[0]>>var[1]>>var[2];
				infile>>var[3]>>var[4]>>var[5];
				infile>>var[6]>>var[7];
				infile>>var[8]>>var[9];
				infile>>var[10];
				infile>>var[11];

				if(infile.eof()) break;				// end of file, break

				if(var[10] > 5) {
					x=var[1]*capa_in; // injected charge in pC
					y=var[8];
					xErr=0;
					yErr=var[9];
				
					g->SetPoint(n,x,y);
					g->SetPointError(n,xErr,yErr);
					n++;
					if(DEBUG) {
						cout << var[1] << endl;
					}
				}

				if(var[0] != last_dac) {
					g->SetMarkerSize(1.0);
					g->SetMarkerStyle(20+style_count);
//					g->SetMarkerColor(2+style_count);
//					g->SetLineColor(2+style_count);
					g->SetMarkerColor(colors[style_count] + style_diff_count);
					g->SetLineColor(colors[style_count] + style_diff_count);

					leg->AddEntry(g,Form("E_thr = %.0f, I_bias = %2.0f",last_ethr,last_dac),"lep");
					//leg->AddEntry(g,Form("E_thr = %.0f, I_bias = %2.0f",160.0-10.0*i,last_dac),"lep");
//					if(i == 0) {
//						low_pass = "with low pass filter";
//					}
//					else {
//						low_pass = "without low pass filter";
//					}
//					leg->AddEntry(g,Form("%s, I_bias = %2.0f",low_pass, last_dac),"lep");
//
//					leg->AddEntry(g,Form("I_bias = %2.0f", last_dac),"lep");

					graph_list->push_back(g);
					g = new TGraphErrors();
					style_count++;
					n=0;
					last_dac = var[0];
					if(var[11] != last_ethr) {
						style_count = 0;
						last_ethr = var[11];
						style_diff_count++;
					}
				}

			}
			// the last set of data
			g->SetMarkerSize(1.0);
			g->SetMarkerStyle(20+style_count);
			g->SetMarkerColor(colors[style_count] + style_diff_count);
			g->SetLineColor(colors[style_count] + style_diff_count);

			leg->AddEntry(g,Form("E_thr = %.0f, I_bias = %2.0f",var[11],var[0]),"lep");
			//leg->AddEntry(g,Form("E_thr = %.0f, I_bias = %2.0f",160.0-10.0*i,var[0]),"lep");
			//leg->AddEntry(g,Form("without low pass filter, I_bias = %2.0f",var[0]),"lep");
			//leg->AddEntry(g,Form("I_bias = %2.0f",var[0]),"lep");

			graph_list->push_back(g);
		}
		infile.close();

	}


	TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
	TMultiGraph *mg = new TMultiGraph();
//	c1->SetGrid();

	nEntr = graph_list->size();
	cout << "Number of plots: " << nEntr << endl;
	for(i=0;i<nEntr;++i)
	{
		mg->Add(graph_list->back());
		graph_list->pop_back();
	}

	if(DEBUG) {
		mg->Print();
	}

	mg->Draw("ap");
//	mg->SetTitle("ToT vs. C");
	mg->GetXaxis()->SetTitle("C_{in}[pC]");
	mg->GetYaxis()->SetTitle("ToT [s]");
//	mg->GetXaxis()->CenterTitle();
//	mg->GetYaxis()->CenterTitle();
	leg->Draw();
	//c1->BuildLegend();
	c1->Update();

	wait();
	delete c1;
	//c1->Close();
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
