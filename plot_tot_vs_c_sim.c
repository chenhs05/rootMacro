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

void plot_tot_vs_c_sim()	{
	const EColor colors[] = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink};
	std::list<char *> * fileName = new std::list<char *>;
	std::list<TGraphErrors *> * graph_list = new std::list<TGraphErrors *>;
	int n,buf;
	double var[16];

	//fileName->push_back("/home/huangshan/gitRepos/ipython_notebook/result/sim_tot_vs_q/20151127_schematics_tot_vs_q_330pF_no_low_pass_filter_val_ibias_val2_ethr_conv.txt");
	fileName->push_back("/home/huangshan/gitRepos/ipython_notebook/result/sim_tot_vs_q/20151128_schematics_tot_vs_q_33pF_no_low_pass_val_ibias_val2_ethr_conv.txt");
	//fileName->push_back("/home/huangshan/gitRepos/ipython_notebook/result/sim_tot_vs_q/20151131_schematics_tot_vs_q_33pF_with_low_pass_val_ibias_val2_ethr_conv.txt");

	ifstream infile;
	int nEntr = fileName->size();

	cout << "##### nEntr of fileName is " << nEntr << endl;

	TLegend* leg = new TLegend(0.1, 0.6, 0.3, 0.9);
	//char * low_pass;
	char * bias_current;

	int style_count=0;
	int style_diff_count = 0;
	for(int i=0;i<nEntr;i++)
	{
		TGraphErrors *g;
		TH1D *h;
		//for the list
		char *last_file_name = fileName->back();
		fileName->pop_back();

		// cout << "fileName: " << last_file_name << endl;
		// cout << "now fileName has nEntr: " << fileName->size() << endl;

		infile.open(last_file_name,ios::in);
		if(!infile) {
			cout<<"\n"<<last_file_name<<" doesn't exist!!"<<"\n"<<"\n";
			return;
		}
		else {
			float x,y,xErr,yErr;
			float x2,y2;
			int style_count=0;
			int style_diff_count = 0;
			float capa_in=33;
			float last_ibias,last_ethr;

			// read in the first line to get the initial value of last_ibias and last_ethr
			infile.seekg(0,ios::beg);
			while(true) {
				infile>>var[0]>>var[1];
				infile>>var[2]>>var[3];
				infile>>var[4]>>var[5]>>var[6];

				if(infile.eof()) break;				// end of file, break
				if(var[0] == 1 ) continue;

				//if(var[0] == 1 ) continue;

				last_ibias=var[0];
				last_ethr=var[1];
				break;
			}
			last_ibias=5;


			n = 0;
			g = new TGraphErrors();

			infile.seekg(0,ios::beg);
			while(true) {
				infile>>var[0]>>var[1];
				infile>>var[2]>>var[3];
				infile>>var[4]>>var[5]>>var[6];

				if(infile.eof()) break;		// end of file, break

				if(var[0] == 1 ) continue;
				if(var[1] != 120) continue;
				if(var[0] != 5 && var[0] != 9 && var[0] != 13) continue;
				if(var[2] > 1.5 ) continue;

				if(1) {
					x=var[2]*capa_in; // injected charge in pC
					y=var[4] * 1e6;
					xErr=0;
					yErr=0;

					if(var[0] != last_ibias) {
						g->SetMarkerSize(0.75);
						g->SetMarkerStyle(20+style_count);
						g->SetMarkerColor(colors[style_count] + style_diff_count);
						g->SetLineColor(colors[style_count] + style_diff_count);

					//	leg->AddEntry(g,Form("E_thr = %.0f, I_bias = %02.0f",last_ethr,last_ibias),"lep");
						if(last_ibias == 5) {
							bias_current = "low discharge current";
						}
						else if(last_ibias == 9) {
							bias_current = "medium discharge current";
						}
						else if(last_ibias == 13) {
							bias_current = "high discharge current";
						}
						leg->AddEntry(g,Form("%s", bias_current),"p");

						graph_list->push_back(g);
						g = new TGraphErrors();
						style_count++;
						n=0;
						last_ibias = var[0];
					}
					if(var[1] != last_ethr) {
						style_count = 0;
						last_ethr = var[1];
						style_diff_count++;
					}

					g->SetPoint(n,x,y);
					g->SetPointError(n,xErr,yErr);
					n++;
					if(DEBUG) {
						cout << var[1] << endl;
					}
				}



			}
			// the last set of data
			g->SetMarkerSize(0.75);
			g->SetMarkerStyle(20+style_count);
			g->SetMarkerColor(colors[style_count] + style_diff_count);
			g->SetLineColor(colors[style_count] + style_diff_count);

		//	leg->AddEntry(g,Form("E_thr = %.0f, I_bias = %02.0f",last_ethr,last_ibias),"lep");
			if(last_ibias == 5) {
				bias_current = "low discharge current";
			}
			else if(last_ibias == 9) {
				bias_current = "medium discharge current";
			}
			else if(last_ibias == 13) {
				bias_current = "high discharge current";
			}
			leg->AddEntry(g,Form("%s", bias_current),"p");

			graph_list->push_back(g);
			style_count ++;
		}
		infile.close();

	}


	TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
	TMultiGraph *mg = new TMultiGraph();

	nEntr = graph_list->size();
	//cout << "Number of plots: " << nEntr << endl;
	for(i=0;i<nEntr;++i)
	{
		mg->Add(graph_list->back());
		graph_list->pop_back();
	}

	if(DEBUG) {
		mg->Print();
	}

	mg->Draw("ap");

	mg->SetTitle("ToT vs. Q_{in}");
	mg->GetXaxis()->SetTitle("Q_{in} [pC]");
	mg->GetYaxis()->SetTitle("ToT [a.u.]");

	leg->Draw();

	c1->Update();


	wait();
	delete c1;
	delete fileName;
	delete graph_list;
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
		if(input){
			done =kTRUE;
		}
	}while(!done);
}
