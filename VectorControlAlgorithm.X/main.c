
// PIC32MZ2048EFH100 Configuration Bit Settings

// DEVCFG3
// USERID = No Setting
#pragma config FMIIEN = ON              // Ethernet RMII/MII Enable (MII Enabled)
#pragma config FETHIO = ON              // Ethernet I/O Pin Select (Default Ethernet I/O)
#pragma config PGL1WAY = ON             // Permission Group Lock One Way Configuration (Allow only one reconfiguration)
#pragma config PMDL1WAY = ON            // Peripheral Module Disable Configuration (Allow only one reconfiguration)
#pragma config IOL1WAY = ON             // Peripheral Pin Select Configuration (Allow only one reconfiguration)
#pragma config FUSBIDIO = ON            // USB USBID Selection (Controlled by the USB Module)

// DEVCFG2
#pragma config FPLLIDIV = DIV_1         // System PLL Input Divider (1x Divider)
#pragma config FPLLRNG  = RANGE_13_26_MHZ// System PLL Input Range (13-26 MHz Input)
#pragma config FPLLICLK = PLL_POSC      // System PLL Input Clock Selection (POSC is input to the System PLL)
#pragma config FPLLMULT = MUL_32     // System PLL Multiplier (PLL Multiply by 24)
#pragma config FPLLODIV = DIV_8         // System PLL Output Clock Divider (2x Divider)
#pragma config UPLLFSEL = FREQ_24MHZ    // USB PLL Input Frequency Selection (USB PLL input is 24 MHz)

// DEVCFG1
#pragma config FNOSC = SPLL             // Oscillator Selection Bits (System PLL)
#pragma config DMTINTV = WIN_127_128    // DMT Count Window Interval (Window/Interval value is 127/128 counter value)
#pragma config FSOSCEN = OFF            // Secondary Oscillator Enable (Disable SOSC)
#pragma config IESO = OFF               // Internal/External Switch Over (Disabled)
#pragma config POSCMOD = EC             // Primary Oscillator Configuration (External clock mode)
#pragma config OSCIOFNC = ON            // CLKO Output Signal Active on the OSCO Pin (Enabled)
#pragma config FCKSM = CSDCMD           // Clock Switching and Monitor Selection (Clock Switch Disabled, FSCM Disabled)
#pragma config WDTPS = PS1048576        // Watchdog Timer Postscaler (1:1048576)
#pragma config WDTSPGM = STOP           // Watchdog Timer Stop During Flash Programming (WDT stops during Flash programming)
#pragma config WINDIS = NORMAL          // Watchdog Timer Window Mode (Watchdog Timer is in non-Window mode)
#pragma config FWDTEN = OFF             // Watchdog Timer Enable (WDT Disabled)
#pragma config FWDTWINSZ = WINSZ_25     // Watchdog Timer Window Size (Window size is 25%)
#pragma config DMTCNT = DMT31           // Deadman Timer Count Selection (2^31 (2147483648))
#pragma config FDMTEN = OFF             // Deadman Timer Enable (Deadman Timer is disabled)

// DEVCFG0
#pragma config DEBUG = OFF              // Background Debugger Enable (Debugger is disabled)
#pragma config JTAGEN = OFF             // JTAG Enable (JTAG Disabled)
#pragma config ICESEL = ICS_PGx1        // ICE/ICD Comm Channel Select (Communicate on PGEC1/PGED1)
#pragma config TRCEN = ON               // Trace Enable (Trace features in the CPU are enabled)
#pragma config BOOTISA = MIPS32         // Boot ISA Selection (Boot code and Exception code is MIPS32)
#pragma config FECCCON = OFF_UNLOCKED   // Dynamic Flash ECC Configuration (ECC and Dynamic ECC are disabled (ECCCON bits are writable))
#pragma config FSLEEP = OFF             // Flash Sleep Mode (Flash is powered down when the device is in Sleep mode)
#pragma config DBGPER = PG_ALL          // Debug Mode CPU Access Permission (Allow CPU access to all permission regions)
#pragma config SMCLR = MCLR_NORM        // Soft Master Clear Enable bit (MCLR pin generates a normal system Reset)
#pragma config SOSCGAIN = GAIN_2X       // Secondary Oscillator Gain Control bits (2x gain setting)
#pragma config SOSCBOOST = ON           // Secondary Oscillator Boost Kick Start Enable bit (Boost the kick start of the oscillator)
#pragma config POSCGAIN = GAIN_2X       // Primary Oscillator Gain Control bits (2x gain setting)
#pragma config POSCBOOST = ON           // Primary Oscillator Boost Kick Start Enable bit (Boost the kick start of the oscillator)
#pragma config EJTAGBEN = NORMAL        // EJTAG Boot (Normal EJTAG functionality)

// DEVCP0
#pragma config CP = OFF                 // Code Protect (Protection Disabled)

// SEQ3

// DEVADC0

// DEVADC1

// DEVADC2

// DEVADC3

// DEVADC4

// DEVADC7

// #pragma config statements should precede project file includes.
// Use project enums instead of #define for ON and OFF.

#include <xc.h>

#include <stdio.h>
#include <stdlib.h>
#include <p32xxxx.h>
#include <stdint.h>
#include <sys/attribs.h>
#include <proc/p32mz2048efh100.h>
#include <math.h>

//FOR 20 MHz EXTERNAL OSC
#define FP 240000000
#define FPBCLK2 80000000
#define CLKPERBAUD 4
#define BAUDRATE 20000000
#define UBRG (FPBCLK2/(CLKPERBAUD*BAUDRATE))-1
#define NRXWORD 32  //change N words RX here
#define NTXWORD 32  //change N words TX here
#define NBYTERX NRXWORD*2+2
#define NBYTETX NTXWORD*2+2
#define NBYTEUSEFUL_RX NRXWORD*2
#define NBYTEUSEFUL_TX NRXWORD*2
#define LED LATGbits.LATG6
#define EN_MK LATGbits.LATG0

#define DT 0.0002
#define PIDEV3 M_PI/3 // pi/3
#define PIMULT2 M_PI*2 //2*pi
#define TWODEVSQRT3 1.15470054 // 2/sqrt(3))
#define SQRT3 sqrt(3)
#define RAD2DEG 180/M_PI
#define RPS2RPM 9.54929659
#define RPM2RPS 1/RPS2RPM

#define KFLTR 256
#define K_RPM_SCALE 1/10
#define K_REF_FLUX_SCALE 1/100
#define K_PI_COEF_SCALE 1/100
#define K_PHASE_SVM_SCALE 16
#define K_UM_SVM_SCALE 8000
#define MAX_SVM_PHASE 360*16
#define K_CURRENT_SCALE 10
#define K_TORQUE_SCALE 10
#define K_FLUX_SCALE 100


#define K_IA_CALIBRATION 1
#define K_IB_CALIBRATION 1
#define K_IC_CALIBRATION 1
#define K_UDC_CALIBRATION 1

#define LR 0.1056
#define RR 0.44
#define LM 0.1
#define ZP 2
#define RPMNOM (3000*RPM2RPS)/ZP
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

//variables for DMA
volatile unsigned char __attribute__((coherent)) DataRx[NBYTERX], DataTx[NBYTETX];
volatile unsigned int __attribute__((coherent)) BlockCRC;
unsigned short int __attribute__((coherent)) DataRx16[NRXWORD], DataTx16[NTXWORD];


uint16_t SystemStates, Ld, ErrCnt;

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



void FpbclkInit (void)
{
    
    uint8_t int_status;
    // Disable global interrupt
    int_status  = __builtin_get_isr_state();
    __builtin_disable_interrupts();
    
    /* Unlock */
    SYSKEY = 0x00000000;
    SYSKEY = 0xAA996655;
    SYSKEY = 0x556699AA;
    
    PB2DIVCLR = 0x8000;
    
    while (PB2DIVbits.PBDIVRDY == 0);
    
    PB2DIVbits.PBDIV = 2; // div = 3 fpbclk2 = 80 Mhz (UART)
    PB2DIVbits.ON = 1;
    
    while (PB3DIVbits.PBDIVRDY == 0);
        
    PB3DIVbits.PBDIV = 47; // div = 48 fpbclk3 = 5 Mhz (Timer))
    PB3DIVbits.ON = 1;
    
    //switch to 240 MHz
    SPLLCONbits.PLLMULT = 24;
    SPLLCONbits.PLLODIV = 2;
    while (!CLKSTATbits.DIVSPLLRDY);
    
    /* Lock */
    SYSKEY = 0x33333333;
     
    __builtin_enable_interrupts();
}

void PinConfig (void) 
{
    ANSELF = 0x00; //F pins as digital
    ANSELG = 0x00; //G pins as digital
    TRISFbits.TRISF0 = 0; // RF0 - as output
    TRISFbits.TRISF1 = 1; // RF1 - as input
    TRISGbits.TRISG6 = 0; //pin led as output
    TRISGbits.TRISG0 = 0; //pin as output. FPGA pin - 49 (Cyclone II)
    U1RXR = 0x0004; // F1 to RX1. FPGA pin - 51 (Cyclone II)
    RPF0R = 0x0001; // F0 to TX1  FPGA pin - 52 (Cyclone II))
}

void UARTInit (void)
{
    U1BRG = UBRG;
    U1MODEbits.BRGH = 1; // high speed x4
     
    //TX
    U1STAbits.UTXISEL = 0b10; // Interrupt is generated and asserted while the transmit buffer contains at least one empty space
    IPC28bits.U1TXIP = 1; // tx interrupt prior 1
    IPC28bits.U1TXIS = 0; // tx interrupt subprior 0
          
    //RX
    U1STAbits.URXISEL = 0b00; // Interrupt flag bit is asserted while receive buffer is not empty 
    IPC28bits.U1RXIP = 2; // rx interrupt prior 2
    IPC28bits.U1RXIS = 0; // rx interrupt subprior 0
    
    IEC3bits.U1EIE = 0;   //disable all error interrupts
    IEC3bits.U1TXIE = 0;  // tx interrupt enable 
    IEC3bits.U1RXIE = 0;  // rx interrupt enable    
    
    U1STAbits.UTXEN = 1; // enable tx
    U1STAbits.URXEN = 1; // enable rx
    
    IFS3bits.U1TXIF = 0; //clear interrupt tx flag
    IFS3bits.U1RXIF = 0; //clear interrupt rx flag
    
    U1MODEbits.ON = 1;  // enable UART module

}

void Timer2Init (void) //timeout timer for DMA
{
    T2CON = 0x00;
    TMR2 = 0;
    PR2 = 0xc350; //10 ms timeout
    
    IPC2bits.T2IP = 6;
    IPC2bits.T2IS = 0;
    IFS0bits.T2IF = 0;
    IEC0bits.T2IE = 1;
    
    T2CONbits.ON = 1;
}

__inline__ unsigned int __attribute__((always_inline)) _VirtToPhys(const void* p) 
{ 
 return (int)p<0?((int)p&0x1fffffffL):(unsigned int)((unsigned char*)p+0x40000000L); 
}


void DMAInit (void)
{
    //DMA CRC INIT
    // ****************
    
    DCRCXOR=0x8005; // Use the Modbus RTU 16 pol 0x8005 - reflected 0xa001
    IEC4bits.DMA2IE = 0; //disable dma ch2 interrupt
    IFS4bits.DMA2IF = 0; //disable dma ch2 interrupt flag
    DCH2CON = 0x00000003; // ch2 priority 3
    DCH2ECON = 0; // no start irqs, no match enabled
    
    
    //TRANSMITION UART INIT
    // ****************
    
    // ch0 - TX, ch1 - RX, ch2 - CRC
    IEC4bits.DMA0IE = 0; //disable dma ch0 interrupt
    IFS4bits.DMA0IF = 0; //clear dma ch0 interrupt flag
    IEC4bits.DMA1IE = 0; //disable dma ch1 interrupt
    IFS4bits.DMA1IF = 0; //disable dma ch1 interrupt flag

    DCH0CON = 0x00000002; // ch0 priority 2 
    DCH1CON = 0x00000002; // chl priority 2
    
    DCH0ECONbits.CHSIRQ = _UART1_TX_VECTOR;
    DCH0ECONbits.SIRQEN = 1;
    
    DCH1ECONbits.CHSIRQ = _UART1_RX_VECTOR;
    DCH1ECONbits.SIRQEN = 1;
    
    //Channel 0
    DCH0SSA = _VirtToPhys(DataTx); // DMA source
    DCH0DSA = _VirtToPhys((void*)&U1TXREG); //DMA destin.
    
    DCH0SSIZ = NBYTETX; //source data block size in bytes
    DCH0DSIZ = 1; 
    DCH0CSIZ = 1;
    
    //Channel 1
    DCH1SSA = _VirtToPhys((void*)&U1RXREG);
    DCH1DSA = _VirtToPhys(DataRx);
        
    DCH1SSIZ = 1;
    DCH1DSIZ = NBYTERX; //destin. data block size in bytes
    DCH1CSIZ = 1;
    
    
    //Interrupts settings
    DCH0INTCLR=0x00ff00ff; //clear all interrupts and events
    DCH0INTbits.CHBCIE = 0; //interrupt event - block transer
       
    IPC33bits.DMA0IP = 0b011; // prior 3
    IPC33bits.DMA0IS = 0b00;  // subprior 0
    IEC4bits.DMA0IE = 0;     // enable interrupt 
    
    DCH1INTCLR=0x00ff00ff;
    DCH1INTbits.CHBCIE = 1;
    
    IPC33bits.DMA1IP = 0b100; // prior 4
    IPC33bits.DMA1IS = 0b00;
    IEC4bits.DMA1IE = 1;
    
    //DCH1CONbits.CHAEN = 1; //auto enable ch1
    DCH1CONbits.CHEN = 1; // enable ch1
    
    DMACONbits.ON = 1; // DMA ON
}

void UARTRXDummy(void)
{
    uint8_t dumm;
    
    while(!U1STAbits.RIDLE);
    U1STAbits.OERR = 0;
    while(U1STAbits.URXDA)
        dumm = U1RXREG;       
}

void TransmitInit (void)
{
    UARTRXDummy();
    DCH0CONbits.CHEN = 1;
    DCH1CONbits.CHEN = 1;
}

unsigned short int Crc16DMA (unsigned char *data, unsigned short int length)
{
    unsigned char PenuByte, LastByte, poll;
    unsigned short int crc;
   
    //save last 2 bytes
    PenuByte = data[length-2];
    LastByte = data[length-1];
    
    //last 2 bytes must be 0x00 0x00 for HW CRC (not counted)
    data[length-2] = 0;
    data[length-1] = 0;
    
    DCRCCON=0x01000fc2; // CRC enabled, polynomial length 16, dma ch - 2, LSb
    DCH2SSA = _VirtToPhys(data); //source
    DCH2DSA = _VirtToPhys(&BlockCRC); 
    DCRCDATA=0xeaa8; // seed the CRC generator 0xeaa8 - reflected 0xFFFF
    DCH2SSIZ = length; //sourse size
    DCH2CSIZ = length; //cell size     
    DCH2DSIZ = 4;
    
    DCH2CONbits.CHEN = 1; //enable ch2 (CRC channel)
    DCH2ECONbits.CFORCE = 1;
    while(!DCH2INTbits.CHBCIF); //successful transfer
    //delay
    poll = length>>1;
    while(poll--); 
    
    //bit-wise swap
    crc = DCRCDATA;
    crc = (((crc & 0xAAAA)>>1) | ((crc & 0x5555)<<1));
    crc = (((crc & 0xCCCC)>>2) | ((crc & 0x3333)<<2));
    crc = (((crc & 0xF0F0)>>4) | ((crc & 0x0F0F)<<4));
    crc = ((crc>>8) | (crc<<8));
    
    //restore last 2 bytes
    data[length-2] = PenuByte;
    data[length-1] = LastByte;
    
    return crc;
}

void ReOrderByteDMA (void *DataSource, void *DataDestination, unsigned char NByte)
{
    unsigned char poll;
    
    DCRCCON = 0x38000f82; // byte reorder
    DCH2SSA = _VirtToPhys(DataSource); //source
    DCH2DSA = _VirtToPhys(DataDestination); //destin
    DCH2SSIZ = NByte; //sourse size
    DCH2CSIZ = NByte; //cell size
    DCH2DSIZ = NByte;
    
    DCH2CONbits.CHEN = 1; //enable ch2 (CRC channel)
    DCH2ECONbits.CFORCE = 1;
    while(!DCH2INTbits.CHBCIF); //successful transfer
    //delay
    poll = NByte>>2;;
    while(poll--);
}


//*** Unsigned Intereger16 MOD ***//

unsigned short UInt16Mod (short X, short Mod)
{
    if(X >= Mod)
        X = X - Mod;
    if(X < -Mod)
        X = X + Mod;
    if(X < 0)
        X = Mod + X;
    
    return ((unsigned short)X);
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
     RefValues.Wr = 0;
}

//*** Set Initial Values ***//

void SetInitValues (void)
{
     MagnetizationPI.Ki = KIMAGNET;
     MagnetizationPI.Kp = KPMAGNET;
     MagnetizationPI.IntegralSaturation = LIMITINTMAGNET;
     MagnetizationPI.OutputSaturation = LIMITOUTMAGNET;
     UdcFltr.K = KFLTR;
     ZeroCalibrationValues.ZeroUdc = INITIALZEROUDC;
}

//*** StartAlgorithm ***//

void StartAlgorithm (void)
{
    ClrVariables();
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
            
            //RefValues.Wr = 0;
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

//*** Data Recieve Conversion ***//

void DataRxConversion (unsigned short *DataRx)
{
    double RefSpeedRPM;
    
    SystemStates = (SystemStates >> 8) << 8 | (char)DataRx[0];
    
    if(SystemStates CHECK REVERS)
        RefSpeedRPM = (double)DataRx[1] *  RPM2RPS * (-1);
    else
        RefSpeedRPM = (double)DataRx[1] *  RPM2RPS;
    RefValues.Wr = SpeedRateLimiter(RefSpeedRPM, RefValues.Wr, DataRx[2], RPMNOM, DT);
    
    if(SystemStates CHECK RUNCONTROL && ~SystemStates CHECK ALGORITHMRUN)
    {
        ZeroCalibrationValues.ZeroIa = DataRx[3];
        ZeroCalibrationValues.ZeroIb = DataRx[4];
    }
    MotorMeasure.Ia = (double)((short)DataRx[3] - ZeroCalibrationValues.ZeroIa) * K_IA_CALIBRATION;
    MotorMeasure.Ib = (double)((short)DataRx[4] - ZeroCalibrationValues.ZeroIb) * K_IB_CALIBRATION;
    MotorMeasure.Ic = -MotorMeasure.Ia - MotorMeasure.Ib;
    
    if(SystemStates CHECK ZEROCALIBRATION)
        ZeroCalibrationValues.ZeroUdc = DataRx[5];
    UdcFltr.In = (double)((short)DataRx[5] - ZeroCalibrationValues.ZeroUdc) * K_UDC_CALIBRATION;
    DigitFltr(&UdcFltr);
    MotorMeasure.Udc = UdcFltr.Out;
    if(MotorMeasure.Udc < 1)
        MotorMeasure.Udc = 1;
    
    MotorMeasure.Wr = (double)((short)DataRx[6]) * RPM2RPS * K_RPM_SCALE;
    RefValues.Flux = (double)DataRx[7] * K_REF_FLUX_SCALE;
    
    SpeedPI.Kp = (double)DataRx[8] * K_PI_COEF_SCALE;
    SpeedPI.Ki = (double)DataRx[9] * K_PI_COEF_SCALE;
    SpeedPI.OutputSaturation = DataRx[10];
    SpeedPI.IntegralSaturation = SpeedPI.OutputSaturation;
    
    FluxPI.Kp = (double)DataRx[11] * K_PI_COEF_SCALE;
    FluxPI.Ki = (double)DataRx[12] * K_PI_COEF_SCALE;
    FluxPI.OutputSaturation = DataRx[13];
    FluxPI.IntegralSaturation = FluxPI.OutputSaturation;
    
    TorquePI.Kp = (double)DataRx[14] * K_PI_COEF_SCALE;
    TorquePI.Ki = (double)DataRx[15] * K_PI_COEF_SCALE;
    TorquePI.OutputSaturation = (double)DataRx[16];
    TorquePI.IntegralSaturation = TorquePI.OutputSaturation;
    
    IdPI.Kp = (double)DataRx[17] * K_PI_COEF_SCALE;
    IdPI.Ki = (double)DataRx[18] * K_PI_COEF_SCALE;
    IdPI.OutputSaturation = DataRx[19];
    IdPI.IntegralSaturation = IdPI.OutputSaturation;
    
    IqPI.Kp = (double)DataRx[20] * K_PI_COEF_SCALE;
    IqPI.Ki = (double)DataRx[21] * K_PI_COEF_SCALE;
    IqPI.OutputSaturation = DataRx[22];
    IqPI.IntegralSaturation = IqPI.OutputSaturation;  
}

//*** Data Transmit Conversion ***//

void DataTxConversion (unsigned short *DataTx)
{
    DataTx[0] = SystemStates >> 8;
    DataTx[1] = ControlValues.Um * K_UM_SVM_SCALE;
    DataTx[2] = UInt16Mod( (ControlValues.Phase * RAD2DEG * K_PHASE_SVM_SCALE), MAX_SVM_PHASE);
    DataTx[3] = MotorMeasure.Ia * K_CURRENT_SCALE;
    DataTx[4] = MotorMeasure.Ib * K_CURRENT_SCALE;
    DataTx[5] = MotorMeasure.Ic * K_CURRENT_SCALE;
    DataTx[6] = MotorMeasure.Udc;
    DataTx[7] = RefValues.Wr * RPS2RPM;
    DataTx[8] = ObserverValues.Torque * K_TORQUE_SCALE;
    DataTx[9] = ObserverValues.PsiR * K_FLUX_SCALE;
    DataTx[10] = ObserverValues.ThetaPsiR * RAD2DEG;
    DataTx[11] = RefValues.Torque * K_TORQUE_SCALE;
    DataTx[12] = RefValues.TorquePI * K_TORQUE_SCALE;
    DataTx[13] = RefValues.PsiRPI * K_FLUX_SCALE;
    DataTx[14] = ControlValues.Usd * K_CURRENT_SCALE;
    DataTx[15] = ControlValues.Usq * K_CURRENT_SCALE;  
}



int main (int argc, char** argv) 
{
    FpbclkInit();
    PinConfig();
    INTCONbits.MVEC = 1; //multivector interrupt mode enable
    UARTInit();
    DMAInit();
    Timer2Init(); //DMA timeout
    EN_MK = 1;
    
    SetInitValues();
    
    while(1)
    {
    
    }
    
    return (EXIT_SUCCESS);
}


void __ISR(_DMA1_VECTOR, ipl4SOFT) _DMA1handler (void) //rx interrupt
{
    uint16_t crc16;
    
    crc16 = Crc16DMA(DataRx, NBYTERX);
    
    if(crc16 == ((DataRx[NBYTERX-2]<<8)|DataRx[NBYTERX-1]))
    {
        TMR2 = 0;
        
        ReOrderByteDMA(DataRx, DataRx16, NBYTEUSEFUL_RX);
        if(Ld >= 1000)
        {
            Ld = 0;
            LED=~LED;
        }
        else
            Ld++;
        
        DataRxConversion(&DataRx16); //read data from DataRx16[n];
        MainAlgorithm(); //execution of the algorithm;
        DataTxConversion(&DataTx16); //write refresh data to DataTx16[n];
    }
    else
    {
        ErrCnt++;   //error processing
    }
    
    DataTx16[NTXWORD-1] = ErrCnt;
    ReOrderByteDMA(DataTx16, DataTx, NBYTEUSEFUL_TX);
    crc16 = Crc16DMA(DataTx, NBYTETX);
    DataTx[NBYTETX-1] = crc16;
    DataTx[NBYTETX-2] = crc16>>8;
  
    T2CONbits.ON = 1;
    DCH1INTCLR = 0x00FF;
    IFS4bits.DMA1IF = 0; 
    TransmitInit();
}

//void __ISR(_DMA0_VECTOR, ipl3SOFT) _DMA0handler (void) //tx interrupt
//{ 
//    DCH0INTCLR=0x00FF;
//    IFS4bits.DMA0IF = 0;
//}

void __ISR(_TIMER_2_VECTOR, ipl6AUTO) _Timer2handler (void) //timeout DMA reset and reenable rx ch1;
{
    uint8_t poll = 0xff;
    
    if(EN_MK)
    {
       EN_MK = 0;
       DCH1CONbits.CHEN = 0;
    }
    else
    {
       DCH1INTCLR = 0x00FF;
       IFS4bits.DMA1IF = 0;
       UARTRXDummy();
       DCH1CONbits.CHEN = 1;
       while(poll--);
       EN_MK = 1;
       T2CONbits.ON = 0;
    }

    TMR2 = 0;
    
    IFS0bits.T2IF = 0;
}