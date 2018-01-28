#include "cpvPitchShift.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "fftw3.h"


#define M_PI 3.14159265358979323846
#define MAX_FRAME_LENGTH 8192

struct cpvStruct{
    double gInFIFO[MAX_FRAME_LENGTH];
    double gFFTworksp[4*MAX_FRAME_LENGTH];
    double gFFTworksp_real[2*MAX_FRAME_LENGTH];
    double interpOutput[2*MAX_FRAME_LENGTH];
    double gLastPhase[MAX_FRAME_LENGTH/2+1];
    double sphase[MAX_FRAME_LENGTH/2+1];
    double gOutputAccum[4*MAX_FRAME_LENGTH];

    double expct;
    long inFifoLatency, fftFrameSize2;

    long stepSize;
    long osamp;
    double sampleRate;
    long fftFrameSize;
    double window[MAX_FRAME_LENGTH];
    fftw_plan pfft;
    fftw_plan pifft;

    bool firsttime;
};

cpvStruct cpvData;

// -----------------------------------------------------------------------------------------------------------------
void interpft(int ny);

void cpvPitchShiftInit(long stepSize_i, long fftFrameSize_i, double sampleRate_i){
    cpvData.firsttime = true;
    
    cpvData.stepSize = stepSize_i;
    cpvData.osamp = fftFrameSize_i/stepSize_i;
    cpvData.sampleRate = sampleRate_i;
    
    cpvData.fftFrameSize = fftFrameSize_i;
	cpvData.fftFrameSize2 = cpvData.fftFrameSize/2;
	cpvData.expct = 2.*M_PI*(double)cpvData.stepSize/(double)cpvData.fftFrameSize;
	cpvData.inFifoLatency = cpvData.fftFrameSize-cpvData.stepSize;
    
    memset(cpvData.gInFIFO, 0, MAX_FRAME_LENGTH*sizeof(double));
    memset(cpvData.gFFTworksp, 0, 2*MAX_FRAME_LENGTH*sizeof(double));
    memset(cpvData.gLastPhase, 0, (MAX_FRAME_LENGTH/2+1)*sizeof(double));
    memset(cpvData.sphase, 0, (MAX_FRAME_LENGTH/2+1)*sizeof(double));
    memset(cpvData.gOutputAccum, 0, 2*MAX_FRAME_LENGTH*sizeof(double));

    for(int k = 0; k < cpvData.fftFrameSize;++k){
        cpvData.window[k] = -.5*cos(2.*M_PI*(double)k/(double)cpvData.fftFrameSize)+.5;
    }

    //fftw_destroy_plan(pfft);
    //fftw_destroy_plan(pfft);
    //pfft = fftw_plan_dft_1d(fftFrameSize, (fftw_complex *)gFFTworksp, (fftw_complex *)gFFTworksp, FFTW_FORWARD, FFTW_ESTIMATE);
    //pifft = fftw_plan_dft_1d(fftFrameSize, (fftw_complex *)gFFTworksp, (fftw_complex *)gFFTworksp, FFTW_BACKWARD, FFTW_ESTIMATE);

    cpvData.pfft = fftw_plan_dft_r2c_1d(cpvData.fftFrameSize, cpvData.gFFTworksp_real, (fftw_complex *)cpvData.gFFTworksp,FFTW_MEASURE);
    cpvData.pifft = fftw_plan_dft_c2r_1d(cpvData.fftFrameSize, (fftw_complex *)cpvData.gFFTworksp, cpvData.gFFTworksp_real,FFTW_MEASURE);
}

void cpvPitchShift(double pitchShift, double *indata, double *outdata)
{
    double magn, phase, tmp, real, imag;
    long k, qpd, index;
    
    if(pitchShift>2) pitchShift=2;
    if(pitchShift<0.5) pitchShift=0.5;
    
    int ifftFrameSize = lround(cpvData.fftFrameSize/pitchShift);
    double Hopratio = cpvData.fftFrameSize/(double)ifftFrameSize;

//     double Hopratio = round(stepSize * pitchShift)/(double)stepSize;
//     int ifftFrameSize = round(fftFrameSize/Hopratio);

    memcpy(&cpvData.gInFIFO[cpvData.inFifoLatency], indata, (size_t) cpvData.stepSize*sizeof(double));
    
    /* do windowing and re,im interleave */
    for (k = 0; k < cpvData.fftFrameSize;k++) {
        cpvData.gFFTworksp_real[k] = cpvData.window[k]*cpvData.gInFIFO[k];
    }

    /* ***************** ANALYSIS ******************* */
    /* do transform */
    fftw_execute(cpvData.pfft);
    
    /* this is the analysis step */
    for (k = 0; k <= cpvData.fftFrameSize2; k++) {

        /* de-interlace FFT buffer */
        real = cpvData.gFFTworksp[2*k];
        imag = cpvData.gFFTworksp[2*k+1];

        /* compute magnitude and phase */
        magn = sqrt(real*real + imag*imag);
        phase = atan2(imag,real);
        
        //magn = 10*log(1+magn/4);

        /* compute phase difference */
        tmp = (phase - cpvData.gLastPhase[k]);
        cpvData.gLastPhase[k] = phase;

        /* subtract expected phase difference */
        tmp -= (double)k*cpvData.expct;
        
        //tmp = tmp - round(tmp/(2*M_PI))*2*M_PI;

        /* map delta phase into +/- Pi interval */
        qpd = tmp/M_PI;
        if (qpd >= 0) qpd += qpd&1;
        else qpd -= qpd&1;
        tmp -= M_PI*(double)qpd;

        
        tmp = (tmp + (double)k*cpvData.expct) * Hopratio;
        
        if(cpvData.firsttime){
            cpvData.sphase[k] = phase;
            cpvData.firsttime =false;
        }else{
            cpvData.sphase[k] += tmp;
        }
        
        cpvData.gFFTworksp[2*k] = magn*cos(cpvData.sphase[k]);
        cpvData.gFFTworksp[2*k+1] = magn*sin(cpvData.sphase[k]);
    }
    
    /* do inverse transform */
    fftw_execute(cpvData.pifft);
    
    for(k=0; k < cpvData.fftFrameSize; k++) {
        cpvData.gFFTworksp_real[k] = cpvData.window[k]*cpvData.gFFTworksp_real[k]/cpvData.fftFrameSize/(cpvData.osamp*3/(double)8)*Hopratio;
    }
    
    interpft(ifftFrameSize);

    if(ifftFrameSize>=cpvData.fftFrameSize){
        int b = (ifftFrameSize-cpvData.fftFrameSize)/2;
        for(k=0; k < cpvData.fftFrameSize; k++) {
            cpvData.gOutputAccum[k] += cpvData.interpOutput[k+b];
        }
    }else{
        int b = (cpvData.fftFrameSize-ifftFrameSize)/2;
        for(k=0; k < ifftFrameSize; k++) {
            cpvData.gOutputAccum[k+b] += cpvData.interpOutput[k];
        }
    }

//     for(k=0; k < ifftFrameSize; k++) {
//         gOutputAccum[k] += interpOutput[k];
//     }
    
    
    for (k = 0; k < cpvData.stepSize; k++) outdata[k] = cpvData.gOutputAccum[k];

    /* shift accumulator */
    memmove(cpvData.gOutputAccum, cpvData.gOutputAccum+cpvData.stepSize, (2*cpvData.fftFrameSize-cpvData.stepSize)*sizeof(double));
    memset(&cpvData.gOutputAccum[2*cpvData.fftFrameSize-cpvData.stepSize-1], 0, cpvData.stepSize*sizeof(double));

    /* move input FIFO */
    for (k = 0; k < cpvData.inFifoLatency; k++) cpvData.gInFIFO[k] = cpvData.gInFIFO[k+cpvData.stepSize];
    
    return;
}

// interpolation function (c translation of matlabs interpft)
void interpft(int ny){
    if(ny==0) return;
    
    int m = cpvData.fftFrameSize;
    
    //If necessary, increase ny by an integer multiple to make ny > m.
    int incr;
    if(ny > m){
        incr = 1;
    }else{
        incr = m/ny + 1;
        ny = incr*ny;
    }
    
    fftw_plan inpifft = fftw_plan_dft_c2r_1d(ny, (fftw_complex *)cpvData.gFFTworksp, cpvData.gFFTworksp_real,FFTW_ESTIMATE);
    
    int nyqst = (m+1)/2+((m+1) % 2 != 0);
    fftw_execute(cpvData.pfft);
    
    memmove(&cpvData.gFFTworksp[2*(nyqst+ny-m)], &cpvData.gFFTworksp[2*nyqst], (m-nyqst)*2*sizeof(double));
    memset(&cpvData.gFFTworksp[2*nyqst], 0, (ny-m)*2*sizeof(double));
    
    if(m%2==0){
        cpvData.gFFTworksp[2*(nyqst-1)] /= 2;
        cpvData.gFFTworksp[2*(nyqst-1)+1] /= 2;
        cpvData.gFFTworksp[2*(nyqst+ny-m-1)] = cpvData.gFFTworksp[2*(nyqst-1)];
        cpvData.gFFTworksp[2*(nyqst+ny-m-1)+1] = cpvData.gFFTworksp[2*(nyqst-1)+1];
    }
    
    fftw_execute(inpifft);
    
    double norm = 1.0/(double) m ; //(double)ny/(m*sqrt((double)m*ny));
    
    for(int i=0;i<ny;++i){
        cpvData.interpOutput[i] = cpvData.gFFTworksp_real[i*incr]*norm;
    }
    
    fftw_destroy_plan(inpifft);
    return;
}

