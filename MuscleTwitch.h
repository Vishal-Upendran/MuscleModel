//Muscle Twitch.h
#ifndef MuscleTwitch_H
#define MuscleTwitch_H
/*Above: Prevent multiple instances of this class*/
class MuscleTwitch
{
	private:
	double *Impulse_response; //Create a lookup table at initialization
	int *spiketrain,window_size; //A spike train of size Delta_T/Sampling time.
	public:
	MuscleTwitch(double*,int); //Constructor
	double GetTwitch(int spike); //Gives the response at current time.
	void UpdateSpike(int); //Updates the spike window
	~MuscleTwitch(){delete[] spiketrain;delete[] Impulse_response;} //Destructor of the system.
};

#endif