//MuscleOutput.h
/* As with all other modules, even thuis has been implemneted as a class. Just inclue the header to start using the module!*/
#ifndef MuscleOutput_H
#define MuscleOutput_H
/*The above two lines to prevent the class from being initialized twice. Chk stackoverflow for more details. */ 
#include "MuscleTwitch.h"
class MuscleOutput: public MuscleTwitch
{
	private:
	double kse,kpe,b,x,xL,x_dot,T; 
	double Fi,Fi_1;
	double srate,Rthresh;
	public:
	MuscleOutput(int,double,double,double,double,double,double []);
	double ForceOutput(double,double,int);	
	void Srate();
};
/*kse: series elastic, kpe: parallel elastic, b: damper, x: muscle length, xL: equilibrium length, x_dot: velocity, T: sampling time, Fi: Force at time t+T, Fi_1: Force at time t, srate: spikerate, Rthresh: threshold spike rate for fiber /
#endif