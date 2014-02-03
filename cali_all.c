gROOT->Reset();
#include <Riostream.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include "math.h"
#include "string.h"
#include "TROOT.h"
#include "Riostream.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"  
#include "TRandom.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPostScript.h"

//Char_t *directory = "E:\\root\\IHep beam test";

//Char_t *directory = "E:\\root\\swap"; //�����ļ��Ĵ��Ŀ¼
//Char_t *fileName="run373"; //�����ļ�����fileName.root

Char_t *directory = "E:\\root\\MTD_data";
Char_t *fileName="035B";

//PMT �� RPC��ʱ�䣬����
Float_t Pm1_t,Pm2_t,Pm3_t,Pm4_t,Pm5_t,Pm6_t;
Float_t q3_1,q3_2,q3_3,q3_4,q1,q2,q3,q4;
Float_t Trpc_t[12],Trpc_q[12];

Char_t buf[1024],buff[1024],psName[1024];

const Int_t ntotmax=500000;
Float_t ttm[ntotmax]; //���������ʱ��
Float_t tot[6][ntotmax];  //��������ķ���
Float_t mean,sigma;
Float_t Nsigma=2.0; //��ʱ��ֲ�����Gaussion fit�ķ�Χ��ȡƽ��ֵ����Nsigma��sigma
const Int_t nCor=5; //����������
const Int_t stripMode=0; //�������ź����ӷ�ʽ��0--0����Ӧ6����1����Ӧ7��... 1--0����Ӧ11����1����Ӧ10��������
Int_t indxHitStrip;  //start from 0�������Ķ��������
Int_t indxHitStrip2;  //���ڴ�Ŷ�������һ�˵����
const Int_t IdxRpcCharge=4;  //RPC�ĵ����tot��ά�����ϵ���ţ���RPC��ɴ�tot[IdxRpcCharge][0]��ʼ��Ĭ��Ϊ4���������Ҫ��������
const Int_t nCharge_pmt=4;
const Int_t nCharge_rpc=5;

/*///IHEP proton
//const Float_t mintot[6]={2800,2400,1500,2800,100,100};     //T-Qͼ��x��ķ�Χ����Ӧ��ͬ�ĵ�ɵķ��ȷ�Χ
//const Float_t maxtot[6]={4100,3900,2600,4100,4100,3200};//T-Qͼ��x��ķ�Χ����Ӧ��ͬ�ĵ�ɵķ��ȷ�Χ
//const Float_t fitmin[6]={2800,2400,1500,2800,100,200};  //�����ķ��ȷ�Χ
//const Float_t fitmax[6]={4100,3900,2600,4100,4100,3000};  //�����ķ��ȷ�Χ
*/
/*///IHEP pion
//const Float_t mintot[6]={500,400,600,300,100,100};     //T-Qͼ��x��ķ�Χ����Ӧ��ͬ�ĵ�ɵķ��ȷ�Χ
//const Float_t maxtot[6]={2600,2200,1500,2400,4100,3200};//T-Qͼ��x��ķ�Χ����Ӧ��ͬ�ĵ�ɵķ��ȷ�Χ
//const Float_t fitmin[6]={500,400,700,300,200,200};  //�����ķ��ȷ�Χ
//const Float_t fitmax[6]={2600,2200,1500,2400,4000,3000};  //�����ķ��ȷ�Χ
*/
/*///��PMTϵͳ
const Float_t mintot[5]={200,200,200,1000,100};     //T-Qͼ��x��ķ�Χ����Ӧ��ͬ�ĵ�ɵķ��ȷ�Χ
const Float_t maxtot[5]={1400,1600,800,2200,4100};		//T-Qͼ��x��ķ�Χ����Ӧ��ͬ�ĵ�ɵķ��ȷ�Χ
const Float_t fitmin[5]={200,200,200,1000,150};  //�����ķ��ȷ�Χ
const Float_t fitmax[5]={1400,1600,800,2200,4100};  //�����ķ��ȷ�Χ
*/
////MTD QC
const Float_t mintot[6]={0,0,200,0,50,200};
const Float_t maxtot[6]={1600,2000,2000,800,4200,3500};
const Float_t fitmin[6]={50,100,250,50,100,500};
const Float_t fitmax[6]={1600,1800,1800,750,4000,3000};

Int_t maxNbins=Int_t((maxtot[4]-mintot[4])/100);  //T-Q��ϵͼx���bin����
Float_t tlimit_pmt=60;    //ʱ��ֲ�ͼ�ķ�Χ�� [-tLimit,tLimit]
Float_t tlimit_rpc=60;    //ʱ��ֲ�ͼ�ķ�Χ�� [-tLimit,tLimit]

/*///well cut - 2 sigma, beam test, proton
//Float_t min_qpair[2]={5550,4584}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4");6463,456.3;5378,396.8
//Float_t max_qpair[2]={7376,6172};
//Float_t min_qdiff[2]={-409,-1977};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4");265.7,337.5;-1409,283.8
//Float_t max_qdiff[2]={951,-841};
//Float_t min_diff[2]={-60,-679}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t");-51.15,4.393;-668.4,5.414
//Float_t max_diff[2]={-42,-658}; 
//Float_t min_position=605; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")618.1,6.741
//Float_t max_position=631;
//Float_t min_tof=425; 	//h1->Draw("(Pm1_t+Pm2_t)-(Pm3_t+Pm4_t)")
//Float_t max_tof=490;
*/
/*/////pion
//Float_t min_qpair[2]={2000,1500}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4");
//Float_t max_qpair[2]={4000,3000};
//Float_t min_qdiff[2]={-300,-1000};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4");
//Float_t max_qdiff[2]={1150,800};
//Float_t min_diff[2]={-70,-690}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t");
//Float_t max_diff[2]={-30,-650}; 
//Float_t min_position=590; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")
//Float_t max_position=650;
//Float_t min_tof=430; 	//h1->Draw("(Pm1_t+Pm2_t)-(Pm3_t+Pm4_t)")
//Float_t max_tof=490;

//Float_t shift_pmt=-137;	//ttm[i]=(Pm1_t-Pm2_t-Pm3_t+Pm4_t)/4+shift_pmt;	
//Float_t shift_rpc=-1633;
*/
/*///run347
//����Cut�ķ�Χ
Float_t min_qpair[2]={900,1500}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4")
Float_t max_qpair[2]={2100,2600};//
Float_t min_qdiff[2]={-700,-1500};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4")
Float_t max_qdiff[2]={500,-700};//
Float_t min_diff[2]={-35,70}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t")
Float_t max_diff[2]={85,200}; //
Float_t min_position=-190; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")
Float_t max_position=10;
//��PMT��RPC��ʱ��ֲ��Ƶ�0������ƫ����
Float_t shift_pmt=184.1; //deltaT=(Pm1_t+Pm2_t-Pm3_t-Pm4_t)/4+shift_pmt
Float_t shift_rpc=-1054.7; //rpc4+rpc10
*/
/*////run373
//����Cut�ķ�Χ
Float_t min_qpair[2]={1000,1200}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4")
Float_t max_qpair[2]={2200,2300};//
Float_t min_qdiff[2]={-700,-1300};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4")
Float_t max_qdiff[2]={500,-300};//
Float_t min_diff[2]={-30,75}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t")
Float_t max_diff[2]={90,200}; //
Float_t min_position=-190; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")
Float_t max_position=-10;
//��PMT��RPC��ʱ��ֲ��Ƶ�0������ƫ����
Float_t shift_pmt=181.9; //deltaT=(Pm1_t+Pm2_t-Pm3_t-Pm4_t)/4+shift_pmt
Float_t shift_rpc=-1095.1; 
*/
/*////023
//����Cut�ķ�Χ
Float_t min_qpair[2]={700,700}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4")
Float_t max_qpair[2]={2300,2000};//
Float_t min_qdiff[2]={-1300,-50};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4")
Float_t max_qdiff[2]={500,1400};//
Float_t min_diff[2]={-85,-25}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t")
Float_t max_diff[2]={60,120}; //
Float_t min_position=-180; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")
Float_t max_position=0;
//��PMT��RPC��ʱ��ֲ��Ƶ�0������ƫ����
Float_t shift_pmt=-41; //deltaT=(Pm1_t+Pm2_t-Pm3_t-Pm4_t)/4+shift_pmt
Float_t shift_rpc[6]={-1160,-1153,-1178,-1145,-1183,-1172}; 
*/
/*//024A
//����Cut�ķ�Χ
Float_t min_qpair[2]={700,600}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4")
Float_t max_qpair[2]={2300,2000};//
Float_t min_qdiff[2]={-1300,-100};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4")
Float_t max_qdiff[2]={600,1300};//
Float_t min_diff[2]={-85,-25}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t")
Float_t max_diff[2]={60,120}; //
Float_t min_position=-180; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")
Float_t max_position=0;
//��PMT��RPC��ʱ��ֲ��Ƶ�0������ƫ����
Float_t shift_pmt=-41; //deltaT=(Pm1_t+Pm2_t-Pm3_t-Pm4_t)/4+shift_pmt
Float_t shift_rpc[6]={-1153,-1067,-1181,-1141,-1166,-1166};
*/ 
/*//024B
//����Cut�ķ�Χ
Float_t min_qpair[2]={700,600}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4")
Float_t max_qpair[2]={2300,2000};//
Float_t min_qdiff[2]={-1300,-100};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4")
Float_t max_qdiff[2]={600,1300};//
Float_t min_diff[2]={-85,-25}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t")
Float_t max_diff[2]={60,120}; //
Float_t min_position=-180; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")
Float_t max_position=0;
//��PMT��RPC��ʱ��ֲ��Ƶ�0������ƫ����
Float_t shift_pmt=-41; //deltaT=(Pm1_t+Pm2_t-Pm3_t-Pm4_t)/4+shift_pmt
Float_t shift_rpc[6]={-1156,-1146,-1178,-1141,-1167,-1167}; 
*/
/*//025A
//����Cut�ķ�Χ
Float_t min_qpair[2]={700,600}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4")
Float_t max_qpair[2]={2300,2000};//
Float_t min_qdiff[2]={-1300,-100};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4")
Float_t max_qdiff[2]={600,1300};//
Float_t min_diff[2]={-85,-25}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t")
Float_t max_diff[2]={60,120}; //
Float_t min_position=-180; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")
Float_t max_position=0;
//��PMT��RPC��ʱ��ֲ��Ƶ�0������ƫ����
Float_t shift_pmt=-41; //deltaT=(Pm1_t+Pm2_t-Pm3_t-Pm4_t)/4+shift_pmt
Float_t shift_rpc[6]={-1159,-1086,-1170,-1146,-1187,-1171}; 
*/
/*//025B
//����Cut�ķ�Χ
Float_t min_qpair[2]={600,600}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4")
Float_t max_qpair[2]={2200,2000};//
Float_t min_qdiff[2]={-1600,-100};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4")
Float_t max_qdiff[2]={300,1300};//
Float_t min_diff[2]={-85,-30}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t")
Float_t max_diff[2]={60,120}; //
Float_t min_position=-180; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")
Float_t max_position=0;
//��PMT��RPC��ʱ��ֲ��Ƶ�0������ƫ����
Float_t shift_pmt=-41; //deltaT=(Pm1_t+Pm2_t-Pm3_t-Pm4_t)/4+shift_pmt
Float_t shift_rpc[6]={-1160,-1072,-1173,-1145,-1185,-1172}; 
*/
/*//020A
//����Cut�ķ�Χ
Float_t min_qpair[2]={800,500}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4")
Float_t max_qpair[2]={2400,1500};//
Float_t min_qdiff[2]={-1100,-100};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4")
Float_t max_qdiff[2]={1100,1200};//
Float_t min_diff[2]={-85,-30}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t")
Float_t max_diff[2]={60,120}; //
Float_t min_position=-200; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")
Float_t max_position=0;
//��PMT��RPC��ʱ��ֲ��Ƶ�0������ƫ����
Float_t shift_pmt=-41; //deltaT=(Pm1_t+Pm2_t-Pm3_t-Pm4_t)/4+shift_pmt
Float_t shift_rpc[6]={-1159,-1072,-1178,-1146,-1186,-1172}; 
*/

/*//020B
//����Cut�ķ�Χ
Float_t min_qpair[2]={600,500}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4")
Float_t max_qpair[2]={1900,1500};//
Float_t min_qdiff[2]={-1300,-0};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4")
Float_t max_qdiff[2]={100,1300};//
Float_t min_diff[2]={-85,-30}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t")
Float_t max_diff[2]={60,120}; //
Float_t min_position=-200; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")
Float_t max_position=0;
//��PMT��RPC��ʱ��ֲ��Ƶ�0������ƫ����
Float_t shift_pmt=-41; //deltaT=(Pm1_t+Pm2_t-Pm3_t-Pm4_t)/4+shift_pmt
Float_t shift_rpc[6]={-1159,-1070,-1175,-1146,-1191,-1172}; 
*/

/*//026A
//����Cut�ķ�Χ
Float_t min_qpair[2]={700,500}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4")
Float_t max_qpair[2]={2000,1500};//
Float_t min_qdiff[2]={-1200,0};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4")
Float_t max_qdiff[2]={500,1300};//
Float_t min_diff[2]={-85,-30}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t")
Float_t max_diff[2]={60,120}; //
Float_t min_position=-200; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")
Float_t max_position=0;
//��PMT��RPC��ʱ��ֲ��Ƶ�0������ƫ����
Float_t shift_pmt=-41; //deltaT=(Pm1_t+Pm2_t-Pm3_t-Pm4_t)/4+shift_pmt
Float_t shift_rpc[6]={-1159,-1157,-1175,-1144,-1182,-1171}; 
*/

/*//026B
//����Cut�ķ�Χ
Float_t min_qpair[2]={600,600}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4")
Float_t max_qpair[2]={2000,1700};//
Float_t min_qdiff[2]={-1200,0};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4")
Float_t max_qdiff[2]={200,1300};//
Float_t min_diff[2]={-85,-30}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t")
Float_t max_diff[2]={60,120}; //
Float_t min_position=-200; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")
Float_t max_position=0;
//��PMT��RPC��ʱ��ֲ��Ƶ�0������ƫ����
Float_t shift_pmt=-41; //deltaT=(Pm1_t+Pm2_t-Pm3_t-Pm4_t)/4+shift_pmt
Float_t shift_rpc[6]={-1159,-1153,-1184,-1145,-1187,-1171}; 
*/

/*//027A
//����Cut�ķ�Χ
Float_t min_qpair[2]={850,600}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4")
Float_t max_qpair[2]={2300,1700};//
Float_t min_qdiff[2]={-1200,0};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4")
Float_t max_qdiff[2]={1200,1400};//
Float_t min_diff[2]={-85,-30}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t")
Float_t max_diff[2]={60,120}; //
Float_t min_position=-200; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")
Float_t max_position=0;
//��PMT��RPC��ʱ��ֲ��Ƶ�0������ƫ����
Float_t shift_pmt=-41; //deltaT=(Pm1_t+Pm2_t-Pm3_t-Pm4_t)/4+shift_pmt
Float_t shift_rpc[6]={-1155,-1148,-1179,-1141,-1170,-1165}; 
*/

/*//035A
//����Cut�ķ�Χ
Float_t min_qpair[2]={750,500}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4")
Float_t max_qpair[2]={2200,1500};//
Float_t min_qdiff[2]={-1200,0};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4")
Float_t max_qdiff[2]={1000,1400};//
Float_t min_diff[2]={-90,-30}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t")
Float_t max_diff[2]={60,120}; //
Float_t min_position=-230; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")
Float_t max_position=20;
//��PMT��RPC��ʱ��ֲ��Ƶ�0������ƫ����
Float_t shift_pmt=-40; //deltaT=(Pm1_t+Pm2_t-Pm3_t-Pm4_t)/4+shift_pmt
Float_t shift_rpc[6]={-1155,-1148,-1179,-1139,-1188,-1167}; 
*/

//035B
//����Cut�ķ�Χ
Float_t min_qpair[2]={450,400}; //h1->Draw("q3_1+q3_2"),h1->Draw("q3_3+q3_4")
Float_t max_qpair[2]={2000,1500};//
Float_t min_qdiff[2]={-1400,0};//h1->Draw("q3_1-q3_2"),h1->Draw("q3_3-q3_4")
Float_t max_qdiff[2]={100,1400};//
Float_t min_diff[2]={-90,-30}; //h1->Draw("Pm1_t-Pm2_t"),h1->Draw("Pm3_t-Pm4_t")
Float_t max_diff[2]={60,100}; //
Float_t min_position=-230; //h1->Draw("(Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)")
Float_t max_position=20;
//��PMT��RPC��ʱ��ֲ��Ƶ�0������ƫ����
Float_t shift_pmt=-40; //deltaT=(Pm1_t+Pm2_t-Pm3_t-Pm4_t)/4+shift_pmt
Float_t shift_rpc[6]={-1155,-1152,-1183,-1144,-1191,-1171}; 


const Float_t totStreamer=1100.0;  //avalanche��streamer�ķֽ��
//for avalanche��avalanche��صĲ���
Float_t mintot_ava[5]={mintot[0],mintot[1],mintot[2],mintot[3],mintot[4]};  //T-Qͼ��x��ķ�Χ����Ӧ��ͬ�ĵ�ɵķ��ȷ�Χ
Float_t maxtot_ava[5]={maxtot[0],maxtot[1],maxtot[2],maxtot[3],totStreamer+100};
Float_t fitmin_ava[5]={fitmin[0],fitmin[1],fitmin[2],fitmin[3],fitmin[4]};		//�����ķ��ȷ�Χ
Float_t fitmax_ava[5]={fitmax[0],fitmax[1],fitmax[2],fitmax[3],totStreamer};
Int_t maxNbins_ava=Int_t((totStreamer+100-mintot[4])/100);  //T-Q��ϵͼx���bin����

//for streamer
Float_t mintot_str[5]={mintot[0],mintot[1],mintot[2],mintot[3],totStreamer-100};		//T-Qͼ��x��ķ�Χ����Ӧ��ͬ�ĵ�ɵķ��ȷ�Χ
Float_t maxtot_str[5]={maxtot[0],maxtot[1],maxtot[2],maxtot[3],maxtot[4]};
Float_t fitmin_str[5]={fitmin[0],fitmin[1],fitmin[2],fitmin[3],totStreamer};		//�����ķ��ȷ�Χ
Float_t fitmax_str[5]={fitmax[0],fitmax[1],fitmax[2],fitmax[3],fitmax[4]};
Int_t maxNbins_str=Int_t((maxtot[4]-totStreamer+100)/100);  //T-Q��ϵͼx���bin����

TF1 *fitf = new TF1("fitf","[0]+[1]/sqrt(x)+[2]/x+[3]*x");  //��������

/******************************************************************************************************************************
//�������� by Huangshan Chen, ��д/��װ�� �������� wjb_test.c
//ttm��ʱ�����ݵ��׵�ַ��һά����Ļ�����������
//tot��������ݵ��׵�ַ����ά���ݵĻ���tot[0];
//nCount: ��Чevent�ĸ����� 
//nCharge����Ҫ�����ĵ�ɸ�����
//nCor������������
//psName�����ps�ļ����ļ����������������directory Ŀ¼������psName_Calibration.ps�ļ���
//tLimit��float��ʽ��ʱ��ֲ�ͼ�ķ�Χ�� [-tLimit,tLimit]��
//tBins: ʱ��ֲ���bin����
//totMin: ���Ȼ���time over threshold�ֲ�����Сֵ�����ڶ�T-Qͼ��x��ķ�Χ����Ӧ��ͬ�ĵ�ɵķ��ȷ�Χ����������׵�ַ��
//totMin: ���Ȼ���time over threshold�ֲ������ֵ�����ڶ�T-Qͼ��x��ķ�Χ����������׵�ַ��
//totBins: T-Q��ϵͼx���bin����
//fitMin�� �����ķ�Χ����Сֵ����������׵�ַ��
//fitMax�� �����ķ�Χ�����ֵ����������׵�ַ��
//fitf�� ����������
//Nsigma�� ��ʱ��ֲ�����Gaussion fit�ķ�Χ��ȡƽ��ֵ����Nsigma��sigma
//nTotMax��ttm����tot[0]�ĳ��ȣ�����ȡ����ȷ��totֵ����Ϊtot[1][0]��tot[0][0]�ĵ�ַ���nTotMax
//IsSplit���Ƿ�RPC�ķ��ȷֳ�avalanche �� streamer �ֱ�����ٷ���һ��fit�õ�ʱ��ֱ��ʣ�0--��1--�ǣ�Ĭ��Ϊ0���������Ҫ������
//totSplit��avalanche��streamer���ȵķֽ�㣬Ĭ��Ϊ0���������Ҫ��������
//IdxRpcCharge��RPC�ĵ����tot��ά�����ϵ���ţ���RPC��ɴ�tot[IdxRpcCharge][0]��ʼ��Ĭ��Ϊ4���������Ҫ��������
******************************************************************************************************************************/
void DoCorrection(Float_t *ttm,Float_t *tot,Int_t nCount,Int_t nCharge,Int_t nCor,Char_t *psName,Float_t tLimit,Float_t tBins,Float_t *totMin,Float_t *totMax,Int_t totBins,Float_t *fitMin,Float_t *fitMax,TF1 *fitf,Float_t Nsigma,Int_t nTotMax,Int_t IsSplit=0,Float_t totSplit=0,Int_t IdxRpcCharge=4)
{
	cout<<"######## "<<psName<<" Calibration ########"<<endl;
	cout<<"Good events:    "<<nCount<<endl;
	
	Double_t mean,sigma,para[8],para2[8],err;
  TH1D *T;
  TH2D *Corr;
  TH1D *htemp;
  Int_t nPad;
  Double_t cor;
  Char_t buf[1024];
  
	TF1 *f1 = new TF1("f1","gaus");
	TF1 *ftemp = new TF1("ftemp","[0]*exp(-0.5*((x-[1])/[2])**2)",-tLimit[0],tLimit[0]);
	f1->SetLineColor(2);
	ftemp->SetLineWidth(3);
	ftemp->SetLineStyle(2);
	ftemp->SetLineColor(2);
	
	sprintf(buf,"%s\\%s_Calibration.ps",directory,psName);
	TPostScript *ps = new TPostScript(buf,112);
  ps->Range(24,16);
  
  TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,800);
	c1->SetFillColor(kWhite);
  c1->GetFrame()->SetFillColor(21);
 	c1->GetFrame()->SetBorderSize(6);
  c1->GetFrame()->SetBorderMode(-1);
  
	gStyle->SetOptFit(110);
	gStyle->SetOptStat(110);
  
  for(Int_t i=0;i<nCor;i++)
  {
  	cout<<"########### The No."<<i+1<<" correction. ##########"<<endl;
  	ps->NewPage();
  	c1->Clear();
  	c1->Divide(1,1);
  	c1->cd(1)->SetLogy();
  	
  	sprintf(buf,"%s_%d_Original",psName,i);
//		sprintf(buf," ");
  	T= new TH1D(buf,buf,tBins,-tLimit,tLimit); 
	  T->SetLineWidth(2);
	  
  	for(Int_t m=0;m<nCount;m++)
  	{
  		T->Fill(ttm[m]);
  	}
  	mean=T->GetMean();
  	sigma=T->GetRMS();
	  f1->SetRange(mean-Nsigma*sigma, mean+Nsigma*sigma);
	  T->Fit("f1","NQR");
	  f1->GetParameters(&para[0]);
	  f1->SetRange(para[1]-Nsigma*para[2], para[1]+Nsigma*para[2]);
	  T->Fit("f1","NQR");
	  f1->GetParameters(&para[0]);
	  f1->SetRange(para[1]-Nsigma*para[2], para[1]+Nsigma*para[2]);
	  T->Draw();
	  T->Fit("f1","QR");
	  f1->GetParameters(&para[0]);
	  
	  ftemp->SetParameters(para[0],para[1],para[2]);
	  ftemp->Draw("same");
	  
//	  T->GetXaxis()->SetTitle("((PMT1+PMT2)-(PMT3+PMT4))/4 [��35ps]");
		T->GetXaxis()->SetTitle("t [35ps]");
		T->GetXaxis()->CenterTitle();
		T->GetYaxis()->SetTitle("Counts");
		T->GetYaxis()->CenterTitle();
		
		c1->cd(1)->Update();
		
  	for(Int_t j=0;j<nCharge;j++)
  	{
  		if(j%2==0) 
  		{
  			ps->NewPage();
  			c1->Clear();
  			c1->Divide(2,2);
  			nPad=1;
  		}
  		else 
  		{
  			nPad=2;
  		}
  		sprintf(buf,"%s_Correction(%d)_%d",psName,i,j);
  		Corr = new TH2D(buf,buf,totBins,totMin[j],totMax[j],tBins,-tLimit,tLimit);
  		
  		for(m=0;m<nCount;m++)
  		{
  			Corr->Fill(*(tot+j*nTotMax+m),ttm[m]);
  		}
  		c1->cd(nPad);
  		Corr->Draw("colz");
  			
  		Corr->GetXaxis()->SetTitle("Charge [0.1pC]");
			Corr->GetXaxis()->CenterTitle();
			Corr->GetYaxis()->SetTitle("t [35ps]");
			Corr->GetYaxis()->CenterTitle();
			c1->cd(nPad)->Update();
  		
  		Corr->ProfileX();
			sprintf(buf,"%s_Correction(%d)_%d_pfx",psName,i,j);
			htemp=(TProfile *)gDirectory->Get(buf);
			c1->cd(nPad)->Update();
			htemp->Draw("QRsame");
			htemp->SetLineColor(kBlack);
			htemp->SetMarkerColor(kBlack);
			htemp->SetMarkerStyle(21);
			htemp->SetMarkerSize(0.6);
  		c1->cd(nPad)->Update();
  		
  		if(IsSplit==1&&j==IdxRpcCharge)
  		{
  			TF1 *fitf2= new TF1("fitf2","fitf");
	  		htemp->Fit("fitf","q","QRsame",fitMin[j],totSplit);
	  		fitf->GetParameters(&para[0]);
	  		fitf2->SetRange(fitMin[j],totSplit);
	  		fitf2->SetParameters(&para[0]);
				fitf2->SetLineWidth(2);
	  		fitf2->Draw("same");
  			c1->cd(nPad)->Update();
	  		htemp->Fit("fitf","q","QRsame",totSplit,fitMax[j]);
	  		fitf->GetParameters(&para2[0]);
  			c1->cd(nPad)->Update();	  		
	  		for(m=0;m<nCount;m++)
	  		{
	  			if((*(tot+j*nTotMax+m))<totSplit)
	  			{
		  			cor=para[0]+para[1]/sqrt((*(tot+j*nTotMax+m)))+para[2]/(*(tot+j*nTotMax+m))+para[3]*(*(tot+j*nTotMax+m));
		  			ttm[m]=ttm[m]-cor;
	  			}
	  			if((*(tot+j*nTotMax+m))>totSplit)
	  			{
		  			cor=para2[0]+para2[1]/sqrt((*(tot+j*nTotMax+m)))+para2[2]/(*(tot+j*nTotMax+m))+para2[3]*(*(tot+j*nTotMax+m));
		  			ttm[m]=ttm[m]-cor;
	  			}
	  		} 				
  		}
  		else
  		{
	  		htemp->Fit("fitf","q","QRsame",fitMin[j],fitMax[j]);
	  		c1->cd(nPad)->Update();
	  		fitf->GetParameters(&para[0]);
	  		for(m=0;m<nCount;m++)
	  		{
	  			cor=para[0]+para[1]/sqrt((*(tot+j*nTotMax+m)))+para[2]/(*(tot+j*nTotMax+m))+para[3]*(*(tot+j*nTotMax+m));
	  			ttm[m]=ttm[m]-cor;
	  		}
  		}
  		
//  		sprintf(buf,"%s-T(%d)_%d",psName,i,j);
  		sprintf(buf," ");
  		T= new TH1D(buf,buf,tBins,-tLimit,tLimit); 
  		for(m=0;m<nCount;m++)
  		{
  			T->Fill(ttm[m]);
  		}
  		c1->cd(nPad+2);
  		c1->cd(nPad+2)->SetLogy();
			c1->cd(nPad+2)->SetBorderSize(6);
	  	mean=T->GetMean();
	  	sigma=T->GetRMS();
		  f1->SetRange(mean-Nsigma*sigma, mean+Nsigma*sigma);
  		T->Fit("f1","NQR");
		  f1->GetParameters(&para[0]);
		  f1->SetRange(para[1]-Nsigma*para[2], para[1]+Nsigma*para[2]);
		  T->Fit("f1","NQR");
		  f1->GetParameters(&para[0]);
		  f1->SetRange(para[1]-Nsigma*para[2], para[1]+Nsigma*para[2]);
		  T->Draw();
		  T->Fit("f1","QR");
		  f1->GetParameters(&para[0]);
		  err=f1->GetParError(2);
		  
		  ftemp->SetParameters(para[0],para[1],para[2]);
		  ftemp->Draw("same");
		  
		  T->GetXaxis()->SetTitle("t [\\times 35ps]");
//		  T->GetXaxis()->SetTitle("((PMT1+PMT2)-(PMT3+PMT4))/4 [\\times 35ps]");
			T->GetXaxis()->CenterTitle();
			T->GetYaxis()->SetTitle("Counts");
			T->GetYaxis()->CenterTitle();
  		c1->cd(nPad+2)->Update();
  	}
  }
  
  ofstream outResult;
  sprintf(buf,"%s\\%s_FitResult.txt",directory,fileName);
  outResult.open(buf,ios::app);
//  outResult<<"########    "<<psName<<"    #######"<<endl;
//  outResult<<"Entries	Sigma	Err of Sigma"<<endl;
  outResult<<indxHitStrip+1<<'\t'<<nCount<<'\t'<<para[2]<<'\t'<<err<<endl;
  outResult.close();
	ps->Close();
	delete c1;
}

//����IHEP����ʵ���ж����������Ƿ�Pion
bool IsPion()
{
	bool isPion;
	if(q3_1>500 && q3_1<2500 && q3_2>500 && q3_2<2000 && q3_3>500 && q3_3<1400 && q3_4>500 && q3_4<2000) isPion=1;
	else isPion=0;
	return isPion;
}
//����IHEP����ʵ���ж����������Ƿ�Proton
bool IsProton()
{
	bool isProton;
	if(q3_1>2900 && q3_1<4000 && q3_2>2500 && q3_2<3800 && q3_3>1600 && q3_3<2500 && q3_4>2900 && q3_4<4000) isProton=1;
	else isProton=0;
	return isProton;
}
//�����ж������Ƿ��������PMT cut֮��
bool InsideCut()
{
	bool isCut=0;
	if(Pm1_t>0 && Pm2_t>0 && Pm3_t>0 &&Pm4_t>0 
		&& q3_1>=10 && q3_2>=10 && q3_3>=10 && q3_4>=10
		&& q3_1+q3_2 >= min_qpair[0] && q3_1+q3_2 <= max_qpair[0]
		&& q3_3+q3_4 >= min_qpair[1] && q3_3+q3_4 <= max_qpair[1]
		&& q3_1-q3_2 >= min_qdiff[0] && q3_1-q3_2 <= max_qdiff[0]
		&& q3_3-q3_4 >= min_qdiff[1] && q3_3-q3_4 <= max_qdiff[1]
		&& Pm1_t-Pm2_t >= min_diff[0] && Pm1_t-Pm2_t <= max_diff[0]
		&& Pm3_t-Pm4_t >= min_diff[1] && Pm3_t-Pm4_t <= max_diff[1]
		&& (Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)>= min_position && (Pm1_t-Pm2_t)-(Pm3_t-Pm4_t)<= max_position
//		&& (Pm1_t+Pm2_t-Pm3_t-Pm4_t)>min_tof && (Pm1_t+Pm2_t-Pm3_t-Pm4_t)<max_tof
		)		{isCut=1;}
		return isCut;
}
//�����ж������Ƿ��Ǻõ��������õ�����ָ�������˶���ʱ���źţ�����֮�ʹ����Ա��������˵ķ���֮�͡��������˶�������������ʵ�������д�ó���
//indxStrip���жϵĶ���������ţ���0��ʼ
//mode���������ź����ӷ�ʽ��0--0����Ӧ6����1����Ӧ7��... 1--0����Ӧ11����1����Ӧ10��������
bool IsGoodEvent(Int_t indxStrip,Int_t mode)
{
	bool isGood=1;
	Int_t indxStrip2;
	if(mode==0)
	{
		if(indxStrip>5) indxStrip=indxStrip-6;
		indxStrip2=indxStrip+6;
	}
	else if(mode==1)
	{
		if(indxStrip>5) indxStrip=11-indxStrip;
		indxStrip2=11-indxStrip;
	}
	
//	Trpc_q[indxStrip]=2*Trpc_q[indxStrip];
	
	if(Trpc_t[indxStrip]==0 || Trpc_t[indxStrip2]==0) 
	{  
		isGood=0;
	}
	else if(mode==0)
	{
		if(indxStrip==0)
		{
			if(Trpc_q[indxStrip]+Trpc_q[indxStrip2]<Trpc_q[indxStrip+1]+Trpc_q[indxStrip2+1])  isGood=0;
		}
		else if(indxStrip==5)
		{
			if(Trpc_q[indxStrip]+Trpc_q[indxStrip2]<Trpc_q[indxStrip-1]+Trpc_q[indxStrip2-1])  isGood=0;
		}
		else if(Trpc_q[indxStrip]+Trpc_q[indxStrip2]<Trpc_q[indxStrip+1]+Trpc_q[indxStrip2+1] || Trpc_q[indxStrip]+Trpc_q[indxStrip2]<Trpc_q[indxStrip-1]+Trpc_q[indxStrip2-1])  isGood=0;
	}
	else if(mode==1)
	{
		if(indxStrip==0)
		{
			if(Trpc_q[indxStrip]+Trpc_q[indxStrip2]<Trpc_q[indxStrip+1]+Trpc_q[indxStrip2-1])  isGood=0;
		}
		else if(indxStrip==5)
		{
			if(Trpc_q[indxStrip]+Trpc_q[indxStrip2]<Trpc_q[indxStrip-1]+Trpc_q[indxStrip2+1])  isGood=0;
		}
		else if(Trpc_q[indxStrip]+Trpc_q[indxStrip2]<Trpc_q[indxStrip+1]+Trpc_q[indxStrip2-1] || Trpc_q[indxStrip]+Trpc_q[indxStrip2]<Trpc_q[indxStrip-1]+Trpc_q[indxStrip2+1])  isGood=0;
	}
	return isGood;
}

//������
void cali_all()
{
	sprintf(buf,"%s\\%s.root",directory,fileName);	//root�ļ���·��
	TFile *file1 = new TFile(buf);
	
	if(file1->IsZombie()) {
    cout<<endl<<buf<<" doesn't exist!!"<<endl<<endl;
    exit(-1);
  }
  else {
  	cout<<"Root File:"<<endl<<"########  "<<buf<<"  ########"<<endl<<endl;
//  	if(stripMode==0)
//		{
//			if(indxHitStrip>5) indxHitStrip=indxHitStrip-6;
//			indxHitStrip2=indxHitStrip+6;
//		}
//		else if(mode==1)
//		{
//			if(indxHitStrip>5) indxHitStrip=11-indxHitStrip;
//			indxHitStrip2=11-indxHitStrip;
//		}
  	ofstream os;
	  sprintf(buf,"%s\\%s_FitResult.txt",directory,fileName);
	  os.open(buf);
	  os<<"########    "<<fileName<<"    ########"<<endl;
	  os<<"indxStrip	Counts	sigma	err of sigma"<<endl;
	  os.close();
	
		sprintf(buf,"h1");  //tree������
		Read(file1,buf);		
		Int_t nEvents= (Int_t)h1->GetEntries();
		
		//all
		DoPmt(0,nEvents);
//		DoRpc(0,nEvents);
		DoRpcSplit(0,nEvents);
		DoRpcAvalanche(0,nEvents);
		DoRpcStreamer(0,nEvents);

//		//first part�����ļ��ֳ��������֣������Ƕ�ǰ�벿�ֵ����ݽ��д���
//		DoPmt(0,Int_t(0.5*nEvents),1);
////		DoRpc(0,Int_t(0.5*nEvents),1);
//		DoRpcSplit(0,Int_t(0.5*nEvents),1);
//		DoRpcAvalanche(0,Int_t(0.5*nEvents),1);
//		DoRpcStreamer(0,Int_t(0.5*nEvents),1);
//
//		//second part,���ļ��ֳ��������֣������ǶԺ�벿�ֵ����ݽ��д���
//		DoPmt(Int_t(0.5*nEvents),nEvents,2);
////		DoRpc(Int_t(0.5*nEvents),nEvents,2);
//		DoRpcSplit(Int_t(0.5*nEvents),nEvents,2);
//		DoRpcAvalanche(Int_t(0.5*nEvents),nEvents,2);
//		DoRpcStreamer(Int_t(0.5*nEvents),nEvents,2);
//
//		for(indxHitStrip=0;indxHitStrip<6;indxHitStrip++)
//  	{
//			if(stripMode==0)
//			{
//				if(indxHitStrip>5) indxHitStrip=indxHitStrip-6;
//				indxHitStrip2=indxHitStrip+6;
//			}
//			else if(mode==1)
//			{
//				if(indxHitStrip>5) indxHitStrip=11-indxHitStrip;
//				indxHitStrip2=11-indxHitStrip;
//			}
//			cout<<"--------  HitStrip:	"<<indxHitStrip+1<<"	"<<indxHitStrip2+1<<"  --------"<<endl;
//			DoRpcSplit(0,nEvents);
//		}
	}
	file1->Close();
}

//����.root�����ļ��������ֱ�ΪTFile���͵ı��������ļ�����tree������
void Read(TFile *f,Char_t *treeName)
{
 	TTree *h1 = f->Get(treeName);
	
	for(Int_t i=0;i<12;i++)
	{
		sprintf(buf,"Trpc_t%d",i+1);
		h1->SetBranchAddress(buf,&Trpc_t[i]);
		sprintf(buf,"Trpc_q%d",i+1);
		h1->SetBranchAddress(buf,&Trpc_q[i]);
	}
   	
   	h1->SetBranchAddress("q3_1",&q3_1);
   	h1->SetBranchAddress("q3_2",&q3_2);
   	h1->SetBranchAddress("q3_3",&q3_3);
   	h1->SetBranchAddress("q3_4",&q3_4);
   	
   	h1->SetBranchAddress("Pm1_t",&Pm1_t);
   	h1->SetBranchAddress("Pm2_t",&Pm2_t);
   	h1->SetBranchAddress("Pm3_t",&Pm3_t);
   	h1->SetBranchAddress("Pm4_t",&Pm4_t);
   	h1->SetBranchAddress("Pm5_t",&Pm5_t);
   	h1->SetBranchAddress("Pm6_t",&Pm6_t);
}
//��PMT���н����������ֱ������ݿ�ʼ��entry��ţ���ֹ��entry��ţ�ģʽ��0--����ȫ��entry��1--����ǰ�벿��entry��2--�����벿entry��
//ע��ģʽ��ѡ��Դ������û��Ӱ�죬ֻ�������γɲ�ͬ��ps�ļ������洦����
void DoPmt(Int_t indxEventDown,Int_t indxEventUp,Int_t mode=0)
{
	Int_t i,j,k;
	Int_t nCount=0;
	Int_t nEvent= (Int_t)h1->GetEntries();
	
	for (i=indxEventDown;i<indxEventUp;i++)
	{
		ttm[i]=0;
		for(Int_t j=0;j<4;j++)
		{
			tot[j][i]=0;
		}
	}
	
	//�����ݣ�������Ч����������
	for (i=indxEventDown;i<indxEventUp;i++)
	{
		h1->GetEntry(i);
//		if(IsPion()==0) continue;
//		if(IsProton()==0) continue;
		if(InsideCut()==0) continue;
		tot[0][nCount]=q3_1;
		tot[1][nCount]=q3_2;
		tot[2][nCount]=q3_3;
		tot[3][nCount]=q3_4;

		ttm[nCount]=(Pm1_t+Pm2_t-Pm3_t-Pm4_t)/4+shift_pmt;
//		ttm[nCount]=(Pm1_t-Pm2_t-Pm3_t+Pm4_t)/4+shift_pmt;	
		nCount++; //counts of valid events
	}
	
	sprintf(psName,"%d_PMT",mode);
	Int_t nCharge=nCharge_pmt;
	Int_t tBins=2*tlimit_pmt;
	DoCorrection(ttm,tot[0],nCount,nCharge,nCor,psName,tlimit_pmt,tBins,mintot,maxtot,maxNbins,fitmin,fitmax,fitf,Nsigma,ntotmax);
}
//��RPC���н�����������avalanche��streamer�źţ������ֱ������ݿ�ʼ��entry��ţ���ֹ��entry��ţ�ģʽ��0--����ȫ��entry��1--����ǰ�벿��entry��2--�����벿entry��
void DoRpc(Int_t indxEventDown,Int_t indxEventUp,Int_t mode=0)
{	
	Int_t i,j,k;
	Int_t nCount;
	Int_t nEvent= (Int_t)h1->GetEntries();
	
	for(Int_t i=0;i<nEvent+2;i++)
	{
		ttm[i]=0;
		for(Int_t j=0;j<IdxRpcCharge+1;j++)
		{
			tot[j][i]=0;
		}
	}

	nCount=0;
	for (i=indxEventDown;i<indxEventUp;i++)
	{
		h1->GetEntry(i);
//		if(IsPion()==0) continue;
//		if(IsProton()==0) continue;
		if(InsideCut()==0) continue;
			
		if(!IsGoodEvent(indxHitStrip,stripMode)) continue;

		tot[4][nCount]=(Trpc_q[indxHitStrip]+Trpc_q[indxHitStrip2])/2.0;
		ttm[nCount]=(Trpc_t[indxHitStrip]+Trpc_t[indxHitStrip2])/2.0-(Pm1_t+Pm2_t+Pm3_t+Pm4_t)/4.0+shift_rpc;		
		tot[0][nCount]=q3_1;
		tot[1][nCount]=q3_2;
		tot[2][nCount]=q3_3;
		tot[3][nCount]=q3_4;

		nCount++; //counts of valid events
	}

	sprintf(psName,"%d_RPC",mode);
	Int_t nCharge=nCharge_rpc;
	Int_t tBins=2*tlimit_rpc;
	DoCorrection(ttm,tot[0],nCount,nCharge,nCor,psName,tlimit_rpc,tBins,mintot,maxtot,maxNbins,fitmin,fitmax,fitf,Nsigma,ntotmax);
}

//��RPC���н�������avalanche��streamer�����ֱ���н����ٷ�һ��fit�õ�ʱ��ֱ��ʣ������ֱ������ݿ�ʼ��entry��ţ���ֹ��entry��ţ�ģʽ��0--����ȫ��entry��1--����ǰ�벿��entry��2--�����벿entry��
void DoRpcSplit(Int_t indxEventDown,Int_t indxEventUp,Int_t mode=0)
{	
	Int_t i,j,k;
	Int_t nCount;
	Int_t nEvent= (Int_t)h1->GetEntries();
	
	for(i=0;i<nEvent+2;i++)
	{
		ttm[i]=0;
		for(Int_t j=0;j<IdxRpcCharge+1;j++)
		{
			tot[j][i]=0;
		}
	}
	
	nCount=0;
	for (i=indxEventDown;i<indxEventUp;i++)
	{
		h1->GetEntry(i);
//		if(IsPion()==0) continue;
//		if(IsProton()==0) continue;
		if(InsideCut()==0) continue;
		if(IsGoodEvent(indxHitStrip,stripMode)==0) continue;

		tot[4][nCount]=(Trpc_q[indxHitStrip]+Trpc_q[indxHitStrip2])/2.0;
//		tot[4][nCount]=(2*Trpc_q[indxHitStrip]+Trpc_q[indxHitStrip2])/2.0;
		ttm[nCount]=(Trpc_t[indxHitStrip]+Trpc_t[indxHitStrip2])/2.0-(Pm1_t+Pm2_t+Pm3_t+Pm4_t)/4.0+shift_rpc[indxHitStrip];
		tot[0][nCount]=q3_1;
		tot[1][nCount]=q3_2;
		tot[2][nCount]=q3_3;
		tot[3][nCount]=q3_4;

		nCount++; //counts of valid events
	}

//	sprintf(psName,"%d_RPC_Split",mode);
	sprintf(psName,"%d_RPC_Split",indxHitStrip+1);
	Int_t nCharge=nCharge_rpc;
	Int_t tBins=2*tlimit_rpc;
	Float_t totSplit=totStreamer;
	Int_t IsSplit=1;
	DoCorrection(ttm,tot[0],nCount,nCharge,nCor,psName,tlimit_rpc,tBins,mintot,maxtot,maxNbins,fitmin,fitmax,fitf,Nsigma,ntotmax,IsSplit,totSplit);
}
//��RPC��avalanche���ֽ��н����������ֱ������ݿ�ʼ��entry��ţ���ֹ��entry��ţ�ģʽ��0--����ȫ��entry��1--����ǰ�벿��entry��2--�����벿entry��
void DoRpcAvalanche(Int_t indxEventDown,Int_t indxEventUp,Int_t mode=0)
{	
	Int_t i,j,k;
	Int_t nCount;
	Int_t nEvent= (Int_t)h1->GetEntries();
	
	for(i=0;i<nEvent+2;i++)
	{
		ttm[i]=0;
		for(Int_t j=0;j<IdxRpcCharge+1;j++)
		{
			tot[j][i]=0;
		}
	}
	
	nCount=0;
	for (i=indxEventDown;i<indxEventUp;i++)
	{
		h1->GetEntry(i);
//		if(IsPion()==0) continue;
//		if(IsProton()==0) continue;
		if(InsideCut()==0) continue;
		if(IsGoodEvent(indxHitStrip,stripMode)==0) continue;
		if(Trpc_q[indxHitStrip]>totStreamer || Trpc_q[indxHitStrip2]>totStreamer) continue;

		tot[4][nCount]=(Trpc_q[indxHitStrip]+Trpc_q[indxHitStrip2])/2.0;
		ttm[nCount]=(Trpc_t[indxHitStrip]+Trpc_t[indxHitStrip2])/2.0-(Pm1_t+Pm2_t+Pm3_t+Pm4_t)/4.0+shift_rpc;		
		tot[0][nCount]=q3_1;
		tot[1][nCount]=q3_2;
		tot[2][nCount]=q3_3;
		tot[3][nCount]=q3_4;

		nCount++; //counts of valid events
	}

	sprintf(psName,"%d_RPC_avalanche",mode);
	Int_t nCharge=nCharge_rpc;
	Int_t tBins=2*tlimit_rpc;
	DoCorrection(ttm,tot[0],nCount,nCharge,nCor,psName,tlimit_rpc,tBins,mintot_ava,maxtot_ava,maxNbins_ava,fitmin_ava,fitmax_ava,fitf,Nsigma,ntotmax);
}
//��RPC��streamer���ֽ��н����������ֱ������ݿ�ʼ��entry��ţ���ֹ��entry��ţ�ģʽ��0--����ȫ��entry��1--����ǰ�벿��entry��2--�����벿entry��
void DoRpcStreamer(Int_t indxEventDown,Int_t indxEventUp,Int_t mode=0)
{	
	Int_t i,j,k;
	Int_t nCount;
	Int_t nEvent= (Int_t)h1->GetEntries();
	
	for (i=indxEventDown;i<indxEventUp;i++)
	{
		ttm[i]=0;
		for(Int_t j=0;j<IdxRpcCharge+1;j++)
		{
			tot[j][i]=0;
		}
	}
	
	nCount=0;
	for (i=indxEventDown;i<indxEventUp;i++)
	{
		h1->GetEntry(i);
//		if(IsPion()==0) continue;
//		if(IsProton()==0) continue;
		if(InsideCut()==0) continue;
		if(IsGoodEvent(indxHitStrip,stripMode)==0) continue;
		if(Trpc_q[indxHitStrip]<totStreamer || Trpc_q[indxHitStrip2]<totStreamer) continue;

		tot[4][nCount]=(Trpc_q[indxHitStrip]+Trpc_q[indxHitStrip2])/2.0;
		ttm[nCount]=(Trpc_t[indxHitStrip]+Trpc_t[indxHitStrip2])/2.0-(Pm1_t+Pm2_t+Pm3_t+Pm4_t)/4.0+shift_rpc;			
		tot[0][nCount]=q3_1;
		tot[1][nCount]=q3_2;
		tot[2][nCount]=q3_3;
		tot[3][nCount]=q3_4;

		nCount++; //counts of valid events
	}

	sprintf(psName,"%d_RPC_streamer",mode);
	Int_t nCharge=nCharge_rpc;
	Int_t tBins=2*tlimit_rpc;
	DoCorrection(ttm,tot[0],nCount,nCharge,nCor,psName,tlimit_rpc,tBins,mintot_str,maxtot_str,maxNbins_str,fitmin_str,fitmax_str,fitf,Nsigma,ntotmax);
}
