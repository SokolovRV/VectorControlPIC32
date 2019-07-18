#include <xc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DT 0.00002
#define PIDEV3 M_PI/3 // pi/3
#define PIMULT2 M_PI*2 //2*pi
#define TWODEVSQRT3 1.15470054 // 2/sqrt(3))

#define RPM2RPS 9.54929659
#define KRPMSCALE 10

#define LR 0.1056
#define RR 0.44
#define LM 0.1
#define ZP 2
#define RPMNOM 3000*KRPMSCALE/ZP
#define LRDEVRR LR/RR
#define RRDECLR RR/LR
#define TORQUECONST 1.5*LM*ZP/LR/RR
#define OMEGACONST RRDECLR * LM

#define INITIALPSIINTEGRAL 0.0001 

#define SET |=
#define CLR &=~
#define CHECK &

#define RUNCONTROL 0x01
#define REVERS 0x02
#define ALGORITHMRUN 0x0100
#define MAGNETIZATION 0x0200

uint16_t SystemStates;

struct PI {
    double Integral;
    double Ki;
    double Kp;
    double IntegralSaturation;
    double OutputSaturation;
    double IntegralInput;
};

struct Measure {
    double Ia;
    double Ib;
    double Ic;
    double Id;
    double Iq;
    double Wr;
} MotorMeasure;

struct Control {
    double Usd;
    double Usq;
    double UsAlpha;
    double UsBeta;
    double Um;
    double Phase;
    double DeltaPsiR;
} ControlValues;

struct Observer {
    double PsiR;
    double PsiRInput;
    double ThetaPsiR;
    double Torque;
    double LimitPsiIntegral;
    double Omega0;
} ObserverValues;

struct Ref {
    double Wr;
    double Flux;
    double PsiRPI;
    double Torque;
    double TorquePI;
    double FluxErr0;
} RefValues;

struct PI SpeedPI;
struct PI FluxPI;
struct PI TorquePI;
struct PI IdPI;
struct PI IqPI;
struct PI MagnetizationPI;


//*** PI regulator ***//

double PIRegulator (double Err, struct PI *PIParam)
{
    double Output;

    PIParam->Integral = PIParam->Integral + PIParam->IntegralInput * DT;   
    if(PIParam->Integral > PIParam->IntegralSaturation)
        PIParam->Integral = PIParam->IntegralSaturation;
    if(PIParam->Integral < -PIParam->IntegralSaturation)
        PIParam->Integral = -PIParam->IntegralSaturation; 
    PIParam->IntegralInput = Err * PIParam->Ki;
    Output = Err * PIParam->Kp + PIParam->Integral;
    if(Output > PIParam->OutputSaturation)
        Output = PIParam->OutputSaturation;
    if(Output < -PIParam->OutputSaturation)
        Output = -PIParam->OutputSaturation;
   
    return (Output);
}


//*** ABC to dq ***//

void ABCtoDq (double ThetaPsiR, struct Measure *MotorMeasure)
{
    MotorMeasure->Id = TWODEVSQRT3 * (MotorMeasure->Ia * sin(ThetaPsiR + PIDEV3) + MotorMeasure->Ib * sin(ThetaPsiR));
    MotorMeasure->Iq = TWODEVSQRT3 * (MotorMeasure->Ia * cos(ThetaPsiR + PIDEV3) + MotorMeasure->Ib * cos(ThetaPsiR));
}

//*** dq to AlphaBeta ***//

void DqtoAlphaBeta (double ThetaPsiR, struct Control *ControlValues)
{
    double CosThetaPsiR, SinThetaPsiR;

    CosThetaPsiR = cos(ThetaPsiR);
    SinThetaPsiR = sin(ThetaPsiR);
    ControlValues->UsAlpha = ControlValues->Usd * CosThetaPsiR - ControlValues->Usq * SinThetaPsiR;
    ControlValues->UsBeta = ControlValues->Usd * SinThetaPsiR + ControlValues->Usq * CosThetaPsiR;
}

//*** Observer / Model Psi R ***//

void Observer (struct Measure *MotorMeasure, struct Observer *ObserverValues)
{
    double WrElectrical;
    
    WrElectrical = MotorMeasure->Wr * ZP;   
    ObserverValues->PsiR = ObserverValues->PsiR + ObserverValues->PsiRInput * DT;    
    if(ObserverValues->PsiR > ObserverValues->LimitPsiIntegral)
        ObserverValues->PsiR = ObserverValues->LimitPsiIntegral;    
    ObserverValues->PsiRInput = (LM * MotorMeasure->Id - ObserverValues->PsiR) * LRDEVRR;      
    ObserverValues->Omega0 = WrElectrical + (MotorMeasure->Iq * OMEGACONST) /  ObserverValues->PsiR ;         
    ObserverValues->Torque = ObserverValues->PsiR * ObserverValues->PsiR * TORQUECONST * (ObserverValues->Omega0 - WrElectrical); 
}


//*** Field-Orient Control ***//

void FOC (void)
{ 
    ObserverValues.ThetaPsiR = ObserverValues.ThetaPsiR + ObserverValues.Omega0 * DT;
    if(ObserverValues.ThetaPsiR >= PIMULT2)
        ObserverValues.ThetaPsiR = ObserverValues.ThetaPsiR - PIMULT2;
    if(ObserverValues.ThetaPsiR <= -PIMULT2)
        ObserverValues.ThetaPsiR = ObserverValues.ThetaPsiR + PIMULT2; 
    ABCtoDq(ObserverValues.ThetaPsiR, &MotorMeasure);
    Observer(&MotorMeasure, &ObserverValues);
    ControlValues.DeltaPsiR = RefValues.Flux - ObserverValues.PsiR;
    RefValues.PsiRPI = PIRegulator(ControlValues.DeltaPsiR, &FluxPI);
    RefValues.Torque = PIRegulator( (RefValues.Wr - MotorMeasure.Wr), &SpeedPI);
    RefValues.TorquePI = PIRegulator( (RefValues.Torque - ObserverValues.Torque), &TorquePI);
    ControlValues.Usq = PIRegulator( (RefValues.TorquePI - MotorMeasure.Iq), &IqPI);
    ControlValues.Usd = PIRegulator( (RefValues.PsiRPI - MotorMeasure.Id), &IdPI);
    DqtoAlphaBeta(ObserverValues.ThetaPsiR, &ControlValues);
    ControlValues.Um = sqrt(ControlValues.UsAlpha * ControlValues.UsAlpha + ControlValues.UsBeta * ControlValues.UsBeta);
    ControlValues.Phase = atan2(ControlValues.UsAlpha, ControlValues.UsBeta);
    if(ControlValues.Um > 1)
        ControlValues.Um = 1;
}

void StartMotor (void)
{
     ObserverValues.PsiR = INITIALPSIINTEGRAL;
     ObserverValues.ThetaPsiR = 0;
     FluxPI.Integral = 0;
     SpeedPI.Integral = 0;
     TorquePI.Integral = 0;
     FluxPI.Integral = 0;
     FluxPI.Integral = 0;
     MagnetizationPI.Integral = 0;
     SystemStates CLR MAGNETIZATION;
     SystemStates SET ALGORITHMRUN;
}

double SpeedRateLimiter (int RefRpm, int TimeNomSpeedInc)
{
    
}

void MainAlgorithm (void)
{
    
}



void main (void)
{
    int a;
    double u, al, be;

    MotorMeasure.Ia = 10.56;
    MotorMeasure.Ib = -5.56;
    MotorMeasure.Wr = 0.1;
    ObserverValues.LimitPsiIntegral = 1000;
    ObserverValues.PsiR = INITIALPSIINTEGRAL;
    RefValues.Flux = 0.8;
    RefValues.Wr = 100;
    
    FluxPI.Ki=2;
    FluxPI.Kp=0.5;
    FluxPI.IntegralSaturation = 10;
    FluxPI.OutputSaturation = 10;
    
    SpeedPI.Ki=2;
    SpeedPI.Kp=0.5;
    SpeedPI.IntegralSaturation = 10;
    SpeedPI.OutputSaturation = 10;
    
    TorquePI.Ki=2;
    TorquePI.Kp=0.5;
    TorquePI.IntegralSaturation = 10;
    TorquePI.OutputSaturation = 10;
    
    IdPI.Ki=2;
    IdPI.Kp=0.5;
    IdPI.IntegralSaturation = 10;
    IdPI.OutputSaturation = 10;
    
    IqPI.Ki=2;
    IqPI.Kp=0.5;
    IqPI.IntegralSaturation = 10;
    IqPI.OutputSaturation = 10;
    
    ObserverValues.ThetaPsiR = -7.2455;
    
    SystemStates = 3;
    if(SystemStates CHECK MAGNETIZATION)
        Nop();
    SystemStates SET MAGNETIZATION;
    Nop();
    while(1)
    {
        Nop();
        for(a=0; a<=1000; a++)

        Nop();
    }
}