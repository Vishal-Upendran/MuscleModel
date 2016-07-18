#include "NuclearChain.h"
/*The way to use is: First, update the model, and then call the primary and secondry firing. */
NuclearChain::NuclearChain(float kp,float ks,float g,float b1,float b2,float t,float a1,float aa)
{
	Gamma=g;kpe=kp;kse=ks;beta1=b1;beta2=b2;T=t;a=a1;alpha=aa;
	xi=xi_1=Fi=Fi_1=0.0;
	xdot=0.0;
}
void NuclearChain::Update(float x,float xd,float g)
{   // If velocity is explicitly given, to reduce computational load.
	xi=x;
	xdot=xd;
	Gamma=g;
	if(xdot>=0)
		b=beta1+beta2*(pow(xdot,a));
	else
		b=beta1-beta2*(pow(pow(xdot,a),(a/2)));
	Fi=((kse*kpe/b)*xi+kse*xdot+kse*Gamma)*T+(1.0-(kpe+kse)*T/b)*Fi_1;
	Fi_1=Fi;
	xi_1=xi;
	Ia=alpha*Fi/kse;
	II=alpha*(xi-Fi/kse);
	
}
void NuclearChain::Update(float x,float g)
{	// If velocity is not given explicitly
	xi=x;
	xdot=(xi-xi_1)/T;
	
	Gamma=g;
	if(xdot>=0)
		b=beta1+beta2*(pow(xdot,a));
	else
		b=beta1-beta2*(pow(pow(xdot,a),(a/2)));
	Fi=((kse*kpe/b)*xi+kse*xdot+kse*Gamma)*T+(1.0-(kpe+kse)*T/b)*Fi_1;
	Fi_1=Fi;
	xi_1=xi;
	Ia=alpha*Fi/kse;
	II=alpha*(xi-Fi/kse);
	
}

/* The functions with Force and all have been defined. Now, we need to add the Ia and II firing rates */
double NuclearChain::Primary()
{
	return Ia;
}
double NuclearChain::Secondary()
{
	return II;
}

