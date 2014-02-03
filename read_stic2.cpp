gROOT->Reset();
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include "math.h"
#include "string.h"
#include "include/myCali.h"
#include "include/EventType.h"
//#include "include/EventDict.h"

char *directory = "/home/huangshan/Documents/measurement/SiPM_Measurements/201310_ctr/";
char *fileName = "dump";

void ReadFile()
{
//	stic_data_t *current_event, *last_event;
//	current_event = new stic_data_t;
//	last_event = new stic_data_t;

	gSystem->Load("include/EventDict.h");
	int chan;
	unsigned short frameNum;

	TFile *fIn = new TFile(Form("%s%s.root",directory,fileName));
	TTree *tree = (TTree*)fIn->Get("dump");
	tree->SetBranchAddress("channel",chan);
	//tree->SetBranchAddress("packet_number",&(current_event->packet_number));
	//tree->SetBranchAddress("frame_number",&(current_event->frame_number));
	//tree->SetBranchAddress("channel",&(current_event->channel));
	//tree->SetBranchAddress("T_CCM",&(current_event->T_CCM));
	//tree->SetBranchAddress("T_CCS",&(current_event->T_CCS));
	//tree->SetBranchAddress("T_badhit",&(current_event->T_badhit));
	//tree->SetBranchAddress("T_fine",&(current_event->T_fine));
	//tree->SetBranchAddress("E_CCM",&(current_event->E_CCM));
	//tree->SetBranchAddress("E_CCS",&(current_event->E_CCS));
	//tree->SetBranchAddress("E_badhit",&(current_event->E_badhit));
	//tree->SetBranchAddress("time",&(current_event->time));
	//tree->SetBranchAddress("energy",&(current_event->energy));
	//tree->SetBranchAddress("errors",&(current_event->errors));
	int nEvents = tree->GetEntries();
	printf("event: %d\n",nEvents);

	nEvents=3;
	for(int i=0;i<nEvents;i++)
	{
		tree->GetEntry(i);
		printf("channel: %d\n",chan);
	}

}


void read_stic2()
{
	char buf[1024];
	sprintf(buf,"%s%s.root",directory,fileName);	
	TFile *f1 = new TFile(buf);
	if(f1->IsZombie())	{
		cout<<"Error openning file!!"<<"\n";
		return;
	}
	else {
		cout<<"File:"<<"\n"<<"########  "<<buf<<"  ########"<<"\n\n";
		f1->Close();
		ReadFile(); 
		//PlotSpectrum();
		//DoCali();
	}
}
