//DynamicBag.h
#ifndef DynamicBag_H
#define DynamicBag_H
//Above line to prevent multiple instances of the class.
/*The function and the structure of the thre spindle fibers is basically the same. Only the parameters differ, and not all have Primary afferents. Non linear summation of signals has NOT been implemented.*/

/*The way to use is: First, update the model, and then call the primary and secondry firing. */
class DynamicBag
{
	float Gamma,kpe,kse,b,beta1,beta2,T,a,alpha;
    float xi,xi_1,xdot;
	float Fi,Fi_1;
	float Ia,II;
	public:
	DynamicBag(float,float,float,float,float,float,float); 
	void Update(float,float,float);
	void Update(float,float);	
	double Primary();
	
};
/*Variables have all been defined in the report. Please refer to it. SIa= Ia - primary firing, SII=II -  secondary firing.*/
#endif