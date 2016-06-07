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

void my1DlinearPlot()	{
	const EColor colors[] = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink};
	std::list<char *> * fileName = new std::list<char *>;
	//std::list<TGraphErrors *> * graph_list = new std::list<TGraphErrors *>;
	//std::list<TGraphErrors *> * graph_list2 = new std::list<TGraphErrors *>;
	//std::list<TGraphErrors *> * graph_list3 = new std::list<TGraphErrors *>;

	fileName->push_back("/home/huangshan/gitRepos/ipython_notebook/result/sim_adc_dnl/calibre_C_CC_mix_signal_600mv_-1m_210m_0p1m_vhd_schematics_conv.txt");
	fileName->push_back("/home/huangshan/gitRepos/ipython_notebook/result/sim_adc_dnl/schematics_mix_signal_600mv_-1m_210mv_0p1m_conv.txt");

	ifstream infile;
	int nEntr = fileName->size();
	
	cout << "##### nEntr of fileName is " << nEntr << endl;

	TLegend* leg = new TLegend(0.1, 0.6, 0.3, 0.9);

	TGraphErrors *g;
	TGraphErrors *g2;
	char *last_file_name;

	int n,buf;
	double var[16];

	int style_count=0;
	int style_diff_count = 0;
	float x,y,xErr,yErr;
	float x2,y2,xErr2,yErr2;

	TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
	c1->Divide(1,2);

	for(int i=0;i<nEntr;i++)
	{

		//for the list
		last_file_name = fileName->back();
		fileName->pop_back();

		infile.open(last_file_name,ios::in);
		if(!infile) {
			cout<<"\n"<<last_file_name<<" doesn't exist!!"<<"\n"<<"\n";
			return;
		}
		else {
			n = 0;
			g = new TGraphErrors();
			h = new TH1D("h","pipeline code density", 150, 50, 200);

			infile.seekg(0,ios::beg);
			while(true)	{
				infile>>var[0]>>var[1]>>var[2];
				infile>>var[3];
				infile>>var[4]>>var[5];

				if(infile.eof()) break;		// end of file, break


				x=var[0];
				y=var[2];
				xErr=0;
				yErr=0;

				g->SetPoint(n,x,y);
				g->SetPointError(n,xErr,yErr);
				n++;
				if(DEBUG) {
					cout << var[1] << endl;
				}

			} // while 

			g->SetMarkerSize(0.5);
			g->SetMarkerStyle(20+style_count);
			g->SetMarkerColor(colors[style_count] + style_diff_count);
			g->SetLineColor(colors[style_count] + style_diff_count);

			//g2->SetMarkerSize(0.5);
			//g2->SetMarkerStyle(20+style_count);
			//g2->SetMarkerColor(colors[style_count] + style_diff_count);
			//g2->SetLineColor(colors[style_count] + style_diff_count);

			if(i == 0 ) leg->AddEntry(g,"schematics","lep");
			if(i == 1 ) leg->AddEntry(g,"C+CC","lep");

			style_count ++;


			// start to draw the plots
			c1->cd(1);
			g->Draw("AP");
			g->LeastSquareLinearFit(n,&)


		} // if-else
		infile.close();

	} // for


	TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
	c1->Divide(1,2);
	c1->cd(1);
	TMultiGraph *mg = new TMultiGraph();
	TMultiGraph *mg2 = new TMultiGraph();
	//TMultiGraph *mg3 = new TMultiGraph();
//	c1->SetGrid();

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
//
//	//mg->SetTitle("ToT vs. Qin");
//	//mg->GetXaxis()->SetTitle("Q_{in}[pC]");
//	//mg->GetYaxis()->SetTitle("ToT [ns]");
	mg->SetTitle("ADC value vs. v_{diff}");
	//mg->GetXaxis()->SetTitle("v_{cm} = 0.6V; v_{diff}[V]");
	//mg->GetYaxis()->SetTitle("ADC pipeline SAR part value");
	mg->GetYaxis()->SetTitle("ADC pipeline value");

	//mg->GetYaxis()->SetDecimals();
	//
	leg->Draw();
	//
	c1->cd(2);
	nEntr = graph_list2->size();
	for(i=0;i<nEntr;++i)
	{
		mg2->Add(graph_list2->back());
		graph_list2->pop_back();
	}
	mg2->Draw("AP");
	//mg2->SetTitle("ADC value vs. v_{diff}");
	mg2->GetXaxis()->SetTitle("v_n = 0.6V; v_p = 0.6V + v_{diff}; v_{diff}[V]");
	//mg2->GetYaxis()->SetTitle("ADC pipeline SAR part value");
	mg2->GetYaxis()->SetTitle("ADC pipeline SAR part value");
	c1->Update();

	wait();
	delete c1;
	delete fileName;
	delete graph_list;
	delete graph_list2;
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
