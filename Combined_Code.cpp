/* Here, I have included 1 muscle fiber, 1 GTO, 1 muscle spindle. The muscle fiber needs spikes for communication, and it forms one module. Spindle and GTO form one module, which emit Ia, Ib and II afferent spike rates. So basically, fibers will input spikes, whereas spindle and GTO will give spike rate.*/
/*One more thing: I have included everything as classes here, not as #includables, to make the code contained in one file.*/
/*EDIT 3: The code compiles. But I its frustrating to enter these many parameter values. I will enter, check the result and update */
#include "stdio.h"
#include <iostream>
#include "math.h"

/*************************** Golgi Tendon Organ **********************************/
class GolgiTendonOrgan
{
    private:
    double K,l0,l1,l2,d0,d1,d2;
    double a,A,b,B,T;
    double F_thresh, F_saturate,sbaseline,sSaturate;
    double Fi_1,Fi_2;
    double srate_now,srate_prev,srate_pprev;
    void discrete_initialize(); // Initialize l0,l1..... from a,b....
    public:
    void GTO_initialize(double,double,double,double,double,double,double,double,double,double); // Constructor
    double current_spikerate(double); //Gives the spike rate for the force input.
};
void GolgiTendonOrgan::GTO_initialize(double sb,double sSat,double K1,double F0,double Fn,double slower,double A1, double faster,double B1,double SampleT)
{	//Initializer.
    K=K1,
    F_thresh=F0; //threshold force
    F_saturate=Fn; //saturation force
    a=slower; //slower exponential
    A=A1; //coefficient of a
    b=faster; //faster exp
    B=B1; //coefficient of b
    T=SampleT; //sampling time
	sbaseline=sb; //baseline firing
	sSaturate=sSat; //saturation firing
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
    Fi=F_saturate;
    else if(Fi<=F_thresh)
    Fi=0.0;
    else    
    {	Fi=Fi-F_thresh;
		srate_now=(K*10*(l0*Fi+l1*Fi_1+l2*Fi_2)-d1*srate_prev-d2*srate_pprev)/d0;
	}
    
    Fi_2=Fi_1;
    Fi_1=Fi;
    
    srate_pprev=srate_prev;
    srate_prev=srate_now;
    return srate_now;
}
/*************** Muscle fiber *************************************/
class MuscleTwitch
{
	private:
	double *Impulse_response; //Create a lookup table at initialization
	public:
	int *spiketrain,window_size; //A spike train of size Delta_T/Sampling time.
	MuscleTwitch(double*,int); //Constructor
	double GetTwitch(int spike); //Gives the response at current time.
	void UpdateSpike(int); //Updates the spike window
	~MuscleTwitch(){delete[] spiketrain;delete[] Impulse_response;} //Destructor of the system.
};
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
// Fiber: 
class MuscleOutput: public MuscleTwitch
{
	private:
	double kse,kpe,b,x,xL,x_dot,T; 
	double Fi,Fi_1;
	double srate,Rthresh;
	public:
	MuscleOutput(double,double,double,double,double,double,double,double []);
	double ForceOutput(double,double,int);	
	void Srate();
};
MuscleOutput::MuscleOutput(double rt,double l,double ks,double kp,double b1,double xL0,double deT,double r[]):MuscleTwitch(r,l)
{  /*Variable initialization. MuscleTwitch is also initialized with this.*/
   kse=ks;kpe=kp;b=b1;xL=xL0;
   T=deT; Fi_1=0;   srate=0;
   Rthresh=rt;
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

/******************** Muscle spindle: called as static fiber ****************************/
class StaticBag
{
	float Gamma,kpe,kse,b,beta1,beta2,T,a,alpha;
    float xi,xi_1,xdot;
	float Fi,Fi_1;
	float Ia,II;
	public:
	void Bag_init(float,float,float,float,float,float,float,float);
	void Update(float,float,float);
	void Update(float,float);	
	double Primary();
	double Secondary();
};
void StaticBag::Bag_init(float kp,float ks,float g,float b1,float b2,float t,float a1,float aa)
{
	Gamma=g;kpe=kp;kse=ks;beta1=b1;beta2=b2;T=t;a=a1;alpha=aa;
	xi=xi_1=Fi=Fi_1=0.0;
	xdot=0.0;
}
void StaticBag::Update(float x,float xd,float g)
{	// If velocity given explicitly, to reduce computational load.
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
void StaticBag::Update(float x,float g)
{	//If velocity not given explicitly
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
double StaticBag::Primary()
{
	return Ia;
}
double StaticBag::Secondary()
{
	return II;
}

/*-----*/
GolgiTendonOrgan gto1;
MuscleOutput *fiber1;
StaticBag spindle1;
/*------------------------- All objects defined. Now to write the main program that will actuate them -------------------*/
void muscle_init(GolgiTendonOrgan &,MuscleOutput *,StaticBag &);
void GetSpikeRates(double &,double &,double &,double,double,double,double);

int main()
{
	muscle_init(gto1,fiber1,spindle1);
	double x,xd,g,F,Fout;
	double SIa;
	double SIb;
	double SII;
	int spike;
	/* Initialized all the modules. I don't know how you would call the function, nevertheless, assume spike is my 'spike' input. spike=1, no spike=0 */
	int i=1;
	while(i)
	{
		std::cin>>spike;
		std::cin>>x; // current displacement.
		std::cin>>xd; //current velocity.
		Fout=fiber1->ForceOutput(x,xd,spike);
		std::cin>>g; //current Gamma
		std::cin>>F; //current force on the GTO.
		/* spindle1.Update(x,xd,g);
		SIa=spindle1.Primary();
		SII=spindle1.Secondary();
		SIb=gto1.current_spikerate(F);  These four lines can be replaced by using:*/ 
		GetSpikeRates(SIa,SII,SIb,x,xd,F,g);
		
		/* Now we have spike rates. These functions and classes can be used directly with some spinnaker functions, in which case 'cin' will be replaced by addSpike() or something. Use the pointers however required. Call the above function to store the spike rates in the pointer.*/
		/*One line here for a condition which detects when to shut down the system: set i to 0. For now, it is taken from the user.*/
		std::cin>>i;
		
	}
	return 0;
	
}
void muscle_init(GolgiTendonOrgan &gto,MuscleOutput *mo,StaticBag &fib)
{
	std::cout<<"Enter the initialization for GTO";
	double *val,*resp;
	val=new double[10];
	for(int i=0;i<10;i++)
		std::cin>>val[i];
	gto.GTO_initialize(val[0],val[1],val[2],val[3],val[4],val[5],val[6], val[7],val[8],val[9]); // GTO initialized.
	delete [] val;
	val= new double[7];
	for(int i=0;i<7;i++)
		std::cin>>val[i];
	int i=val[1];
	resp=new double[i];
	for(int i=0;i<val[1];i++)
		std::cin>>resp[i];
		
	fiber1=new MuscleOutput(val[0],val[1],val[2],val[3],val[4],val[5],val[6],resp); //Muscle output initialized.
	delete [] val;
	delete [] resp;
	val=new double[8];
	for(int i=0;i<8;i++)
		std::cin>>val[i];
	fib.Bag_init(val[0],val[1],val[2],val[3],val[4],val[5],val[6],val[7]); // Spindle initialized.
	
	/* The most crudest way of initialization (according to me) has been implemented. */
}
void GetSpikeRates(double &SIa, double &SII, double &SIb, double x, double xd, double F, double g)
{
	spindle1.Update(x,xd,g);
	SIa=spindle1.Primary();
	SII=spindle1.Secondary();
	SIb=gto1.current_spikerate(F);
}
