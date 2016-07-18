//GolgiTendonOrgan.h
#ifndef GolgiTendonOrgan_H
#define GolgiTendonOrgan_H
// Above line to prevent multiple instances of class.
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
    GolgiTendonOrgan(double,double,double,double,double,double,double,double,double); // Constructor
    double current_spikerate(double); //Gives the spike rate for the force input.
};
/* The variables have all been defined in the report. Same names have been used. F_thresh : Threshold force, F_saturate : saturation force. sbaseline: minimum firing rate of GTO. Fi_1, Fi_2: prev 2 times forces, srate_now,etc: Similarly prev and now spike rates.*/
#endif