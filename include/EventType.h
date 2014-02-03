#ifndef EVENT_TYPE_H__
#define EVENT_TYPE_H__

#include "vector"

#include "TObject.h"

class stic_data_t: public TObject
{
public:
	unsigned short			packet_number;	//The Packet ID the event was transmitted in
	unsigned short			frame_number;	//The frame number of the event (for debug)
	unsigned short			channel;	//The channel number
	unsigned int			T_CCM;		//The Time Coarse Counter Master value
	unsigned int			T_CCS;		//The Time Coarse Counter Slave value
	bool				T_badhit;	//The bad hit flag of the time measurement
	unsigned short			T_fine;		//The fine counter value of the time measurement
	unsigned int			E_CCM;		//The coarse counter value of the energy measurement (MASTER)
	unsigned int			E_CCS;		//The coarse counter value of the energy measurement (SLAVE)
	bool				E_badhit;	//The bad hit flag of the energy measurement (not needed here)

	/* higher level entries */
	unsigned int			time;		//time including energy
	unsigned int			energy;		//energy (coarse counter only)			both select the correct CM/CS values
	std::vector<unsigned int>	errors;		//errors reconstruction
ClassDef(stic_data_t,1);
	//void stic_data_t::stic_data_t() {}
};
//ClassImp(stic_data_t);

#endif

