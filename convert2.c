gROOT->Reset();
#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

Char_t gPATH[]="data.txt";
Float_t Trpc_q1,Trpc_q2,Trpc_q3,Trpc_q4,Trpc_q5,Trpc_q6,Trpc_q7,Trpc_q8,Trpc_q9,Trpc_q10,Trpc_q11,Trpc_q12,Trpc_t1,Trpc_t2,Trpc_t3,Trpc_t4,Trpc_t5,Trpc_t6,Trpc_t7,Trpc_t8,Trpc_t9,Trpc_t10,Trpc_t11,Trpc_t12,Pm1_t,Pm2_t,Pm3_t,Pm4_t,Pm5_t,Pm6_t,q3_1,q3_2,q3_3,q3_4;

void convert2()
{
	TFile *f = new TFile("Trpc.root","RECREATE");
	TTree *h1 = new TTree("h1","data from File");
		
	h1->Branch("Trpc_q1",&Trpc_q1,"Trpc_q1/F");
	h1->Branch("Trpc_q2",&Trpc_q2,"Trpc_q2/F");
	h1->Branch("Trpc_q3",&Trpc_q3,"Trpc_q3/F");
	h1->Branch("Trpc_q4",&Trpc_q4,"Trpc_q4/F");
	h1->Branch("Trpc_q5",&Trpc_q5,"Trpc_q5/F");
	h1->Branch("Trpc_q6",&Trpc_q6,"Trpc_q6/F");
	h1->Branch("Trpc_q7",&Trpc_q7,"Trpc_q7/F");
	h1->Branch("Trpc_q8",&Trpc_q8,"Trpc_q8/F");
	h1->Branch("Trpc_q9",&Trpc_q9,"Trpc_q9/F");
	h1->Branch("Trpc_q10",&Trpc_q10,"Trpc_q10/F");
	h1->Branch("Trpc_q11",&Trpc_q11,"Trpc_q11/F");
	h1->Branch("Trpc_q12",&Trpc_q12,"Trpc_q12/F");

	h1->Branch("Trpc_t1",&Trpc_t1,"Trpc_t1/F");
	h1->Branch("Trpc_t2",&Trpc_t2,"Trpc_t2/F");
	h1->Branch("Trpc_t3",&Trpc_t3,"Trpc_t3/F");
	h1->Branch("Trpc_t4",&Trpc_t4,"Trpc_t4/F");
	h1->Branch("Trpc_t5",&Trpc_t5,"Trpc_t5/F");
	h1->Branch("Trpc_t6",&Trpc_t6,"Trpc_t6/F");
	h1->Branch("Trpc_t7",&Trpc_t7,"Trpc_t7/F");
	h1->Branch("Trpc_t8",&Trpc_t8,"Trpc_t8/F");
	h1->Branch("Trpc_t9",&Trpc_t9,"Trpc_t9/F");
	h1->Branch("Trpc_t10",&Trpc_t10,"Trpc_t10/F");
	h1->Branch("Trpc_t11",&Trpc_t11,"Trpc_t11/F");
	h1->Branch("Trpc_t12",&Trpc_t12,"Trpc_t12/F");

	h1->Branch("Pm1_t",&Pm1_t,"Pm1_t/F");
	h1->Branch("Pm2_t",&Pm2_t,"Pm2_t/F");
	h1->Branch("Pm3_t",&Pm3_t,"Pm3_t/F");
	h1->Branch("Pm4_t",&Pm4_t,"Pm4_t/F");
	h1->Branch("Pm5_t",&Pm5_t,"Pm5_t/F");
	h1->Branch("Pm6_t",&Pm6_t,"Pm6_t/F");
	
	h1->Branch("q3_1",&q3_1,"q3_1/F");
	h1->Branch("q3_2",&q3_2,"q3_2/F");
	h1->Branch("q3_3",&q3_3,"q3_3/F");
	h1->Branch("q3_4",&q3_4,"q3_4/F");
	
	ifstream is;
	is.open(gPATH,ios::in);
	while(is>>Trpc_q1>>Trpc_q2>>Trpc_q3>>Trpc_q4>>Trpc_q5>>Trpc_q6>>Trpc_q7>>Trpc_q8>>Trpc_q9>>Trpc_q10>>Trpc_q11>>Trpc_q12>>q3_1>>q3_2>>q3_3>>q3_4>>Trpc_t1>>Trpc_t2>>Trpc_t3>>Trpc_t4>>Trpc_t5>>Trpc_t6>>Trpc_t7>>Trpc_t8>>Trpc_t9>>Trpc_t10>>Trpc_t11>>Trpc_t12>>Pm1_t>>Pm2_t>>Pm3_t>>Pm4_t>>Pm5_t>>Pm6_t)
	{
		if(Pm1_t>0 && Pm2_t>0&& Pm3_t>0&& Pm4_t>0)
		h1->Fill();
//		cout<<Trpc_q1;
	}
	h1->Write();
	is.close();
	f->Close();
	}

