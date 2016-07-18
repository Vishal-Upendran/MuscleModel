//MuscleTwitch.cpp
#include "MuscleTwitch.h"
/*See MuscleTwitch.h for initialization of class and explanation of variables. */

MuscleTwitch::MuscleTwitch(double Response[], int length)
{ /* Constructor function: Initialize impulse response to the given data. And spike train to zeros(1,Delta_T) */
   Impulse_response=new double[length];
   spiketrain=new int[length];
   window_size=length;
   for(int i=0;i<length;i++)
   {
	   Impulse_response[i]=Response[i]; //To be determined by the user. See the sample twitch given in the dataset for an example.
	   spiketrain[i]=0;
   }
}
void MuscleTwitch::UpdateSpike(int spike)
{	//Update the spike train required for convolution. Window size will depend on how fast the impulse response falls. 
	for(int i=0;i<window_size-1;i++)
		spiketrain[i]=spiketrain[i+1]; //pushing value through a queue.
	spiketrain[window_size-1]=spike;
}
double MuscleTwitch::GetTwitch(int spike)
{ // Get the twitch response NOW.
	double Current_response=0.0;
	MuscleTwitch::UpdateSpike(spike); //push the spike onto queue
	for(int i=0;i<window_size;i++)
		Current_response+=spiketrain[i]*Impulse_response[window_size-i]; //Convolution. Take a set of spikes, and convolve with the response. 
	return Current_response;
}