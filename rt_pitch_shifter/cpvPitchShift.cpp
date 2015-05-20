
#include "cpvPitchShift.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "fftw3.h"


#define M_PI 3.14159265358979323846
#define MAX_FRAME_LENGTH 8192


float gInFIFO[MAX_FRAME_LENGTH];
float gFFTworksp[4*MAX_FRAME_LENGTH];
float gFFTworksp_real[2*MAX_FRAME_LENGTH];
float interpOutput[2*MAX_FRAME_LENGTH];
float gLastPhase[MAX_FRAME_LENGTH/2+1];
float sphase[MAX_FRAME_LENGTH/2+1];
float gOutputAccum[4*MAX_FRAME_LENGTH];

double freqPerBin, expct;
long inFifoLatency, fftFrameSize2;

long stepSize;
long osamp;
float sampleRate;
long fftFrameSize;
float window[MAX_FRAME_LENGTH];
fftwf_plan pfft;
fftwf_plan pifft;

bool firsttime = false;

// -----------------------------------------------------------------------------------------------------------------
void interpft(int ny);

void cpvPitchShiftInit(long stepSize_i, long osamp_i, float sampleRate_i){
    firsttime = true;
    
    stepSize = stepSize_i;
    osamp = osamp_i;
    sampleRate = sampleRate_i;
    
    fftFrameSize = stepSize*osamp;
	fftFrameSize2 = fftFrameSize/2;
	freqPerBin = sampleRate/(double)fftFrameSize;
	expct = 2.*M_PI*(double)stepSize/(double)fftFrameSize;
	inFifoLatency = fftFrameSize-stepSize;
    
    memset(gInFIFO, 0, MAX_FRAME_LENGTH*sizeof(float));
    memset(gFFTworksp, 0, 2*MAX_FRAME_LENGTH*sizeof(float));
    memset(gLastPhase, 0, (MAX_FRAME_LENGTH/2+1)*sizeof(float));
    memset(sphase, 0, (MAX_FRAME_LENGTH/2+1)*sizeof(float));
    memset(gOutputAccum, 0, 2*MAX_FRAME_LENGTH*sizeof(float));

    for(int k = 0; k < fftFrameSize;++k){
        window[k] = -.5*cos(2.*M_PI*(double)k/(double)fftFrameSize)+.5;
    }

    //fftwf_destroy_plan(pfft);
    //fftwf_destroy_plan(pfft);
    //pfft = fftwf_plan_dft_1d(fftFrameSize, (fftwf_complex *)gFFTworksp, (fftwf_complex *)gFFTworksp, FFTW_FORWARD, FFTW_ESTIMATE);
    //pifft = fftwf_plan_dft_1d(fftFrameSize, (fftwf_complex *)gFFTworksp, (fftwf_complex *)gFFTworksp, FFTW_BACKWARD, FFTW_ESTIMATE);

    pfft = fftwf_plan_dft_r2c_1d(fftFrameSize, gFFTworksp_real, (fftwf_complex *)gFFTworksp,FFTW_MEASURE);
    pifft = fftwf_plan_dft_c2r_1d(fftFrameSize, (fftwf_complex *)gFFTworksp, gFFTworksp_real,FFTW_MEASURE);
}

void cpvPitchShift(float pitchShift, float *indata, float *outdata)
{
    double magn, phase, tmp, real, imag;
    long k, qpd, index;
    
    if(pitchShift>2) pitchShift=2;
    if(pitchShift<0.5) pitchShift=0.5;
    
    int ifftFrameSize = lround(fftFrameSize/pitchShift);
    double Hopratio = fftFrameSize/(double)ifftFrameSize;

    
    memcpy(&gInFIFO[inFifoLatency], indata, (size_t) stepSize*sizeof(float));
    
    /* do windowing and re,im interleave */
    for (k = 0; k < fftFrameSize;k++) {
        gFFTworksp_real[k] = window[k]*gInFIFO[k];
    }

    /* ***************** ANALYSIS ******************* */
    /* do transform */
    fftwf_execute(pfft);
    
    /* this is the analysis step */
    for (k = 0; k <= fftFrameSize2; k++) {

        /* de-interlace FFT buffer */
        real = gFFTworksp[2*k];
        imag = gFFTworksp[2*k+1];

        /* compute magnitude and phase */
        magn = sqrt(real*real + imag*imag);
        phase = atan2(imag,real);

        /* compute phase difference */
        tmp = (phase - gLastPhase[k]);
        gLastPhase[k] = phase;

        /* subtract expected phase difference */
        tmp -= (double)k*expct;
        
        //tmp = tmp - round(tmp/(2*M_PI))*2*M_PI;

        /* map delta phase into +/- Pi interval */
        qpd = tmp/M_PI;
        if (qpd >= 0) qpd += qpd&1;
        else qpd -= qpd&1;
        tmp -= M_PI*(double)qpd;

        
        tmp = (tmp + (double)k*expct) * Hopratio;
        
        if(firsttime){
            sphase[k] = phase;
            firsttime =false;
        }else{
            sphase[k] += tmp;
        }
        
        gFFTworksp[2*k] = magn*cos(sphase[k]);
        gFFTworksp[2*k+1] = magn*sin(sphase[k]);
    }
    
    /* do inverse transform */
    fftwf_execute(pifft);
    
    for(k=0; k < fftFrameSize; k++) {
        gFFTworksp_real[k] = window[k]*gFFTworksp_real[k]/fftFrameSize/(osamp*3/(double)8);
    }
    
    interpft(ifftFrameSize);

    if(ifftFrameSize>=fftFrameSize){
        int b = (ifftFrameSize-fftFrameSize)/2;
        for(k=0; k < fftFrameSize; k++) {
            gOutputAccum[k] += interpOutput[k+b];
        }
    }else{
        int b = (fftFrameSize-ifftFrameSize)/2;
        for(k=0; k < ifftFrameSize; k++) {
            gOutputAccum[k+b] += interpOutput[k];
        }
    }

//     for(k=0; k < ifftFrameSize; k++) {
//         gOutputAccum[k] += interpOutput[k];
//     }
    
    
    for (k = 0; k < stepSize; k++) outdata[k] = gOutputAccum[k];

    /* shift accumulator */
    memmove(gOutputAccum, gOutputAccum+stepSize, (2*fftFrameSize-stepSize)*sizeof(float));
    memset(&gOutputAccum[2*fftFrameSize-stepSize-1], 0, stepSize*sizeof(float));

    /* move input FIFO */
    for (k = 0; k < inFifoLatency; k++) gInFIFO[k] = gInFIFO[k+stepSize];
    
    return;
}

void interpft(int ny){
    if(ny==0) return;
    
    int m = fftFrameSize;
    
    //If necessary, increase ny by an integer multiple to make ny > m.
    int incr;
    if(ny > m){
        incr = 1;
    }else{
        incr = m/ny + 1;
        ny = incr*ny;
    }
    
    fftwf_plan inpifft = fftwf_plan_dft_c2r_1d(ny, (fftwf_complex *)gFFTworksp, gFFTworksp_real,FFTW_ESTIMATE);
    
    int nyqst = (m+1)/2+((m+1) % 2 != 0);
    fftwf_execute(pfft);
    
    memmove(&gFFTworksp[2*(nyqst+ny-m)], &gFFTworksp[2*nyqst], (m-nyqst)*2*sizeof(float));
    memset(&gFFTworksp[2*nyqst], 0, (ny-m)*2*sizeof(float));
    
    if(m%2==0){
        gFFTworksp[2*(nyqst-1)] /= 2;
        gFFTworksp[2*(nyqst-1)+1] /= 2;
        gFFTworksp[2*(nyqst+ny-m-1)] = gFFTworksp[2*(nyqst-1)];
        gFFTworksp[2*(nyqst+ny-m-1)+1] = gFFTworksp[2*(nyqst-1)+1];
    }
    
    fftwf_execute(inpifft);
    
    float norm = (double)ny/(m*sqrt((double)m*ny));
    
    for(int i=0;i<ny;++i){
        interpOutput[i] = gFFTworksp_real[i*incr]*norm;
    }
    
    return;
}

