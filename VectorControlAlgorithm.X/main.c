#include <xc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DT 1
#define PIDEV3 1.04719755
#define TWODEVSQRT3 1.15470054



struct PI {
    double Integral;
    double Ki;
    double Kp;
    double IntegralSaturation;
    double OutputSaturation;
};

struct Measure {
    double Ia;
    double Ib;
    double Ic;
    double Id;
    double Iq;
    double Wr;
    double Ibrake;
    double Udc;
} MotorMeasure;

struct PI SpeedPI;
struct PI FluxPI;
struct PI TorquePI;
struct PI IdPI;
struct PI IqPI;


//*** PI regulator ***//

double PIRegulator (double Err, struct PI *PIParam)
{
    double Output;
    
    PIParam->Integral = PIParam->Integral + Err * PIParam->Ki * DT;
    
    if(PIParam->Integral > PIParam->IntegralSaturation)
        PIParam->Integral = PIParam->IntegralSaturation;
    if(PIParam->Integral < -PIParam->IntegralSaturation)
        PIParam->Integral = -PIParam->IntegralSaturation;
    
    Output = Err * PIParam->Kp + PIParam->Integral;
       
    if(Output > PIParam->OutputSaturation)
        Output = PIParam->OutputSaturation;
    if(Output < -PIParam->OutputSaturation)
        Output = -PIParam->OutputSaturation;
    
    return (Output);
}


//*** ABC to dq ***//

void ABCtoDq (double Psi, struct Measure *MotorMeasure)
{
    MotorMeasure->Id = TWODEVSQRT3 * (MotorMeasure->Ia * sin(Psi + PIDEV3) + MotorMeasure->Ib * sin(Psi));
    MotorMeasure->Iq = TWODEVSQRT3 * (MotorMeasure->Ia * cos(Psi + PIDEV3) + MotorMeasure->Ib * cos(Psi));
}


void main (void)
{
    double a;

    MotorMeasure.Ia = 5;
    MotorMeasure.Ib = 10;
    MotorMeasure.Ic = -15;
    
    while(1)
    {
        ABCtoDq(2, &MotorMeasure);
        Nop();
    }
}