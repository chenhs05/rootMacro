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
