//Nuclear chain.h
#ifndef NuclearChain_H
#define NuclearChain_H
/*The way to use is: First, update the model, and then call the primary and secondry firing. */
class NuclearChain
{
	float Gamma,kpe,kse,b,beta1,beta2,T,a,alpha;
    float xi,xi_1,xdot;
	float Fi,Fi_1;
	float Ia,II;
	public:
	NuclearChain(float,float,float,float,float,float,float);
	void Update(float,float,float);
	void Update(float,float);	
	double Primary();
	double Secondary();
};
/*Parameters defined in the report. Please go through.*/
#endif