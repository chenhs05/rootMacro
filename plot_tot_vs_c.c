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

void plot_tot_vs_c()	{
	const EColor colors[] = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink};
	std::list<char *> * fileName = new std::list<char *>;
	std::list<TGraphErrors *> * graph_list = new std::list<TGraphErrors *>;
	int n,buf;
	double var[16];

//	fileName->push_back("/home/huangshan/gitRepos/ipython_notebook/result/tot_vs_q/20151119_t_vs_tot_C330p_w_low_pass_filter_conv.txt");
	fileName->push_back("/home/huangshan/gitRepos/ipython_notebook/result/tot_vs_q/20151120_tot_vs_q_C330p_no_low_pass_filter_conv.txt");

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
			float capa_in=330;
			float last_ibias,last_ethr;

			// read in the first line to get the initial value of last_ibias and last_ethr
			infile.seekg(0,ios::beg);
			while(true) {
				infile>>var[0]>>var[1];
				infile>>var[2]>>var[3]>>var[4];
				infile>>var[5]>>var[6]>>var[7];
				infile>>var[8]>>var[9]>>var[10];
				infile>>var[11];

				if(infile.eof()) break;				// end of file, break

				//if(var[0] == 1 ) continue;

				last_ibias=var[0];
				last_ibias=4;
				last_ethr=var[11];
				break;
			}



			n = 0;
			g = new TGraphErrors();

			infile.seekg(0,ios::beg);
			while(true) {
				infile>>var[0]>>var[1];
				infile>>var[2]>>var[3]>>var[4];
				infile>>var[5]>>var[6]>>var[7];
				infile>>var[8]>>var[9]>>var[10];
				infile>>var[11];

				if(infile.eof()) break;		// end of file, break

				if(var[11] == 160) continue;	//omit the ethr==160 cases, too low E_threshold
				if(var[11] == 155) continue;	//omit the ethr==155 cases, too low E_threshold
				if(var[11] != 140) continue;
				if(var[0] != 4 && var[0] != 7 && var[0] != 10) continue;
				//if(var[0] == 1 ) continue;
				if(var[1] > 1.5 ) continue;

				if(var[10] > 5) {
					x=var[1]*capa_in; // injected charge in pC
					y=var[8] * 1e6;
					xErr=0;
					yErr=var[9];
					//x=var[2] * capa_in;
					//y=var[4] * 1e6;

				//	x=var[0];
				//	y=var[2];
				//	x2=var[0];
				//	y2=var[1];
				//	xErr=0;
				//	yErr=0;

				//	if(var[0] >= 0.0003 && var[0] <=0.2062) {
				//		h -> Fill(var[2]);
				//	}

					g->SetPoint(n,x,y);
					g->SetPointError(n,xErr,yErr);
					n++;
					if(DEBUG) {
						cout << var[1] << endl;
					}
				}

				if(var[0] != last_ibias) {
					g->SetMarkerSize(0.75);
					g->SetMarkerStyle(20+style_count);
					g->SetMarkerColor(colors[style_count] + style_diff_count);
					g->SetLineColor(colors[style_count] + style_diff_count);

					//leg->AddEntry(g,Form("I_bias = %02.0f",last_ibias),"lep");
					//leg->AddEntry(g,Form("E_thr = %.0f, I_bias = %02.0f",last_ethr,last_ibias),"lep");
					//if(i == 0) {
					//	low_pass = "with low pass filter";
					//}
					//else {
					//	low_pass = "without low pass filter";
					//}
					//leg->AddEntry(g,Form("%s, I_bias = %2.0f",low_pass, last_ibias),"lep");

					if(last_ibias == 4) {
						bias_current = "low discharge current";
					}
					else if(last_ibias == 7) {
						bias_current = "medium discharge current";
					}
					else if(last_ibias == 10) {
						bias_current = "high discharge current";
					}
					leg->AddEntry(g,Form("%s", bias_current),"lep");

					graph_list->push_back(g);
					g = new TGraphErrors();
					style_count++;
					n=0;
					last_ibias = var[0];

					//if( style_count > 11) {
					//	style_count = 0;
					//	style_diff_count++;
					//}
					if(var[11] != last_ethr) {
						style_count = 0;
						last_ethr = var[11];
						style_diff_count++;
					}
				}

			}
			// the last set of data
			g->SetMarkerSize(0.75);
			g->SetMarkerStyle(20+style_count);
			g->SetMarkerColor(colors[style_count] + style_diff_count);
			g->SetLineColor(colors[style_count] + style_diff_count);


			if(last_ibias == 4) {
				bias_current = "low discharge current";
			}
			else if(last_ibias == 7) {
				bias_current = "medium discharge current";
			}
			else if(last_ibias == 10) {
				bias_current = "high discharge current";
			}
			leg->AddEntry(g,Form("%s", bias_current),"lep");

			//leg->AddEntry(g,Form("I_bias = %2.0f",last_ibias),"lep");
			//leg->AddEntry(g,Form("E_thr = %.0f, I_bias = %2.0f",last_ethr,last_ibias),"lep");
			//leg->AddEntry(g,Form("E_thr = %.0f, I_bias = %2.0f",160.0-10.0*i,var[0]),"lep");
			//
			//if(i == 0 ) leg->AddEntry(g,"schematics","lep");
			//if(i == 1 ) leg->AddEntry(g,"C+CC","lep");

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

	mg->SetTitle("ToT vs. Qin");
	mg->GetXaxis()->SetTitle("Q_{in}[pC]");
	mg->GetYaxis()->SetTitle("ToT [us]");

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
