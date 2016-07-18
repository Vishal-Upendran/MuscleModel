// MuscleOutput.cpp

#include "MuscleOutput.h"
/* Variables details in MuscleOutput.h */
MuscleOutput::MuscleOutput(int rt,int l,double ks,double kp,double b1,double xL0,double deT,double r[]):MuscleTwitch(r,l)
{  /*Variable initialization. MuscleTwitch is also initialized with this.*/
   kse=ks;kpe=kp;b=b1;xL=xL0;
   T=deT; Fi_1=0;   srate=0;
   Rthresh=rt
}
double MuscleOutput::ForceOutput(double x1,double xd,int spike)
{   /*Gives the force output. For the equation, please refer to the report. */
	double Twitch_t=GetTwitch(spike)/Rthresh; //Current twitch
	x=x1;x_dot=xd;
	if(srate>=Rthresh) //Minimum rate required to produce the force.
		Fi=((kse*kpe/b)*(x-xL)+kse*x_dot)*T+(kse/b)*Twitch_t+Fi_1*(1.0-(kse+kpe)*T/b);
	else
		Fi=0;
	Fi_1=Fi;
	return Fi;
}
void MuscleOutput::Srate()
{	// Calculate the spike rate. This is not requried if the Force function is modified to input spike rate separately.
	if(spiketrain[window_size-1]==1)
		srate=(srate*T+spiketrain[window_size-1]-spiketrain[0])/T;
	
}