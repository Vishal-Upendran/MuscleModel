//GolgiTendonOrgan.cpp
#include "GolgiTendonOrgan.h"
/*Check the report and GolgiTendonOrgan.h for variable documentation.*/
GolgiTendonOrgan::GolgiTendonOrgan(double sb,double sSat,double K1,double F0,double Fn,double slower,double A1, double faster,double B1,double SampleT)
{	//Initializer.
    K=K1,
    F_thresh=F0;
    F_saturate=Fn;
    a=slower;
    A=A1;
    b=faster;
    B=B1;
    T=SampleT;
	sbaseline=sb;
	sSaturate=sSat;
    Fi_1=Fi_2=srate_prev=srate_pprev=0.0;
    discrete_initialize();
}
void GolgiTendonOrgan::discrete_initialize()
{  //Refer the report for formulae.
    l0=4.0*(1.0+A+B)+2.0*T*(a+b+A*b+a*B)+a*b*T*T;
    l1=2.0*(a*b*(T*T)-4.0*(1.0+A+B));
    l2=4.0*(1+A+B)-2.0*T*(a+b+A*b+a*B)+a*b*(T*T);
    d0=4.0+T*(2.0*a+2.0*b+a*b*T);
    d1=2.0*a*b*(T*T)-8.0;
    d2=T*(a*b*T-2.0*a-2.0*b)+4.0;
   // cout<<l0<<"\n"<<l1<<"\n"<<l2<<"\n"<<d0<<"\n"<<d1<<"\n"<<d2<<"\n";
}
double GolgiTendonOrgan::current_spikerate(double Fi)
{	/*First step: Saturation and Threshold forces. Force less than threshold causes no change in baseline firing, and force above saturation causes no change in firing rate.*/
    if(Fi>=F_saturate)
    srate_now=sSaturate; 
    else if(Fi<=F_thresh)
    srate_now=sbaseline;
    else    
    srate_now=(K*10*(l0*Fi+l1*Fi_1+l2*Fi_2)-d1*srate_prev-d2*srate_pprev)/d0;
    
    Fi_2=Fi_1;
    Fi_1=Fi;
    
    srate_pprev=srate_prev;
    srate_prev=srate_now;
    return srate_now;
}