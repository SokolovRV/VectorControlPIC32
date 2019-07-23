#include <xc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DT 0.00002
#define PIDEV3 M_PI/3 // pi/3
#define PIMULT2 M_PI*2 //2*pi
#define TWODEVSQRT3 1.15470054 // 2/sqrt(3))
#define SQRT3 sqrt(3)

#define KFLTR 256
#define RPM2RPS 1/9.54929659
#define KRPMSCALE 1/10
#define KFLUXSCALE 1/100
#define K_PI_COEF_SCALE 1/100

#define K_IA_CALIBRATION 1
#define K_IB_CALIBRATION 1
#define K_IC_CALIBRATION 1
#define K_UDC_CALIBRATION 1

#define LR 0.1056
#define RR 0.44
#define LM 0.1
#define ZP 2
#define RPMNOM 3000/ZP
#define LRDEVRR LR/RR
#define RRDECLR RR/LR
#define TORQUECONST 1.5*LM*ZP/LR/RR
#define OMEGACONST RRDECLR * LM
#define FLUXERR0 0.01
#define CURRENTNOM 20

#define KPMAGNET 0.04
#define KIMAGNET 1
#define LIMITINTMAGNET 0.5
#define LIMITOUTMAGNET 0.5

#define INITIALPSIINTEGRAL 0.0001
#define INITIALZEROUDC 2048

#define SET |=
#define CLR &=~
#define CHECK &

#define RUNCONTROL 0x01
#define REVERS 0x02
#define ZEROCALIBRATION 0x04
#define ALGORITHMRUN 0x0100
#define MAGNETIZATION 0x0200

uint16_t SystemStates;

struct Calibration {
    short ZeroIa;
    short ZeroIb;
    short ZeroIc;
    short ZeroUdc;
} ZeroCalibrationValues;

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
    double Udc;
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
} RefValues;

struct DigitFltr {
    double In;
    double Out;
    double K;
    double Accum;
};

struct PI SpeedPI;
struct PI FluxPI;
struct PI TorquePI;
struct PI IdPI;
struct PI IqPI;
struct PI MagnetizationPI;
struct DigitFltr UdcFltr;


//*** Data Recieve Conversion ***//

void DataRxConversion (unsigned short *DataRx)
{
    double RefSpeedRPM;
    
    SystemStates = (SystemStates >> 8) << 8 | (char)DataRx[0];
    
    if(SystemStates CHECK REVERS)
        RefSpeedRPM = (double)DataRx[1] *  RPM2RPS * (-1);
    else
        RefSpeedRPM = (double)DataRx[1] *  RPM2RPS;
    RefValues.Wr = (RefSpeedRPM, RefValues.Wr, DataRx[2], RPMNOM, DT);
    
    if(SystemStates CHECK RUNCONTROL && ~SystemStates CHECK ALGORITHMRUN)
    {
        ZeroCalibrationValues.ZeroIa = DataRx[3];
        ZeroCalibrationValues.ZeroIb = DataRx[4];
    }
    MotorMeasure.Ia = (double)((short)DataRx[3] - ZeroCalibrationValues.ZeroIa) * K_IA_CALIBRATION;
    MotorMeasure.Ib = (double)((short)DataRx[4] - ZeroCalibrationValues.ZeroIa) * K_IB_CALIBRATION;
    MotorMeasure.Ic = -MotorMeasure.Ia - MotorMeasure.Ib;
    
    if(SystemStates CHECK ZEROCALIBRATION)
        ZeroCalibrationValues.ZeroUdc = DataRx[5];
    UdcFltr.In = (double)((short)DataRx[5] - ZeroCalibrationValues.ZeroUdc) * K_UDC_CALIBRATION;
    DigitFltr(&UdcFltr);
    MotorMeasure.Udc = UdcFltr.Out;
    if(MotorMeasure.Udc < 1)
        MotorMeasure.Udc = 1;
    
    MotorMeasure.Wr = (double)DataRx[6] * RPM2RPS * KRPMSCALE;
    RefValues.Flux = (double)DataRx[7] * KFLUXSCALE;
    
    SpeedPI.Kp = (double)DataRx[8] * K_PI_COEF_SCALE;
    SpeedPI.Ki = (double)DataRx[9] * K_PI_COEF_SCALE;
    SpeedPI.OutputSaturation = (double)DataRx[10] * K_PI_COEF_SCALE;
    SpeedPI.IntegralSaturation = SpeedPI.OutputSaturation;
    
    FluxPI.Kp = (double)DataRx[11] * K_PI_COEF_SCALE;
    FluxPI.Ki = (double)DataRx[12] * K_PI_COEF_SCALE;
    FluxPI.OutputSaturation = (double)DataRx[13] * K_PI_COEF_SCALE;
    FluxPI.IntegralSaturation = FluxPI.OutputSaturation;
    
    TorquePI.Kp = (double)DataRx[14] * K_PI_COEF_SCALE;
    TorquePI.Ki = (double)DataRx[15] * K_PI_COEF_SCALE;
    TorquePI.OutputSaturation = (double)DataRx[16] * K_PI_COEF_SCALE;
    TorquePI.IntegralSaturation = TorquePI.OutputSaturation;
    
    IdPI.Kp = (double)DataRx[17] * K_PI_COEF_SCALE;
    IdPI.Ki = (double)DataRx[18] * K_PI_COEF_SCALE;
    IdPI.OutputSaturation = (double)DataRx[19] * K_PI_COEF_SCALE;
    IdPI.IntegralSaturation = IdPI.OutputSaturation;
    
    IqPI.Kp = (double)DataRx[20] * K_PI_COEF_SCALE;
    IqPI.Ki = (double)DataRx[21] * K_PI_COEF_SCALE;
    IqPI.OutputSaturation = (double)DataRx[22] * K_PI_COEF_SCALE;
    IqPI.IntegralSaturation = IqPI.OutputSaturation;  
}


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
    ControlValues.Um = sqrt(ControlValues.UsAlpha * ControlValues.UsAlpha + ControlValues.UsBeta * ControlValues.UsBeta) / MotorMeasure.Udc * SQRT3; ;
    ControlValues.Phase = atan2(ControlValues.UsAlpha, ControlValues.UsBeta);
    if(ControlValues.Um > 1)
        ControlValues.Um = 1;
}

//*** Clear Variables ***//

void ClrVariables (void)
{
     ObserverValues.PsiR = INITIALPSIINTEGRAL;
     ObserverValues.ThetaPsiR = 0;
     FluxPI.Integral = 0;
     SpeedPI.Integral = 0;
     TorquePI.Integral = 0;
     IdPI.Integral = 0;
     IqPI.Integral = 0;
     MagnetizationPI.Integral = 0;
}

//*** Set Initial Values ***//

void SetInitValues (void)
{
     MagnetizationPI.Ki = KIMAGNET;
     MagnetizationPI.Kp = KPMAGNET;
     MagnetizationPI.IntegralSaturation = LIMITINTMAGNET;
     MagnetizationPI.OutputSaturation = LIMITOUTMAGNET;
     SpeedFltr.K = KFLTR;
     ZeroCalibrationValues.ZeroUdc = INITIALZEROUDC;
}

//*** StartAlgorithm ***//

void StartAlgorithm (void)
{
    ClrVariables();
    SetInitValues();
    SystemStates CLR MAGNETIZATION;
    SystemStates SET ALGORITHMRUN;
}

//*** Speed Rate Limiter ***//

double SpeedRateLimiter (double RefSpeed, double ActualRefSpeed, double TimeNomSpeedInc, double NomSpeed, double dT)
{
    double Increment, Delta;
    
    Increment = NomSpeed / (TimeNomSpeedInc / dT);
    Delta = fabs(RefSpeed - ActualRefSpeed);
    if(Delta > Increment)
    {
        if(ActualRefSpeed < RefSpeed)
            ActualRefSpeed = ActualRefSpeed + Increment;
        if(ActualRefSpeed > RefSpeed)
            ActualRefSpeed = ActualRefSpeed - Increment;
    }
    else
        ActualRefSpeed = RefSpeed;
    
    return(ActualRefSpeed);
}

//*** Digital Filtr ***//

void DigitFltr (struct DigitFltr *FltrData)
{
    FltrData->Accum += FltrData->In - FltrData->Out;
    FltrData->Out = FltrData->Accum / FltrData->K;
}

//*** MAIN ALGORITHM ***//

void MainAlgorithm (void)
{
    if(SystemStates CHECK RUNCONTROL && ~SystemStates CHECK ALGORITHMRUN)
        StartAlgorithm();
    if(~SystemStates CHECK RUNCONTROL)
        SystemStates CLR ALGORITHMRUN;
    if(SystemStates CHECK ALGORITHMRUN)
    {
        if(~SystemStates CHECK MAGNETIZATION)
        {
            double CurrentErr;
            
            RefValues.Wr = 0;
            FOC();
            CurrentErr = CURRENTNOM - (sqrt(MotorMeasure.Ia*MotorMeasure.Ia + MotorMeasure.Ib*MotorMeasure.Ib + MotorMeasure.Ic*MotorMeasure.Ic));
            ControlValues.Um = PIRegulator(CurrentErr, &MagnetizationPI);
            ControlValues.Phase = 0;
            if(fabs(ControlValues.DeltaPsiR) <= FLUXERR0)
            {
               FluxPI.Integral = 0;
               SpeedPI.Integral = 0;
               TorquePI.Integral = 0;
               IdPI.Integral = 0;
               IqPI.Integral = 0;
               SystemStates SET MAGNETIZATION;
            }
        }
        else
        {
            FOC();
        }
    }
}



void main (void)
{
    int a;
    double u, al, be;

    MotorMeasure.Ia = 5;
    MotorMeasure.Ib = 5;
    MotorMeasure.Ic = -10;
    MotorMeasure.Wr = 0.1;
    ObserverValues.LimitPsiIntegral = 1000;
    ObserverValues.PsiR = INITIALPSIINTEGRAL;
    RefValues.Flux = 0.8;
    RefValues.Wr = 100;
    MotorMeasure.Udc = 100;
    
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
    
    unsigned short ia=1900;
    short za=2000;
    Nop();
    MotorMeasure.Ia = (double)((short)ia - za) * K_IA_CALIBRATION;
    Nop();
    
    while(1)
    {
        SystemStates SET RUNCONTROL;
        
        for(a=0; a<=1000; a++)
        {
            Nop();
            MainAlgorithm();
            if(a==3)
             RefValues.Flux = 0;    
            Nop();
        }
        Nop();
    }
}