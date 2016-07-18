//Static bag.h

#ifndef StaticBag_H
#define StaticBag_H
/*The way to use is: First, update the model, and then call the primary and secondry firing. */
class StaticBag
{
	float Gamma,kpe,kse,b,beta1,beta2,T,a,alpha;
    float xi,xi_1,xdot;
	float Fi,Fi_1;
	float Ia,II;
	public:
	StaticBag(float,float,float,float,float,float,float);
	void Update(float,float,float);
	void Update(float,float);	
	double Primary();
	double Secondary();
};

#endif