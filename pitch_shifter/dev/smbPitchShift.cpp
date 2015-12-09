#include "smbPitchShift.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "fftw3.h"


#define M_PI 3.14159265358979323846
#define MAX_FRAME_LENGTH 8192

struct smbStruct{
    double gInFIFO[MAX_FRAME_LENGTH];
    double gFFTworksp[2*MAX_FRAME_LENGTH];
    double gFFTworksp_real[MAX_FRAME_LENGTH];
    double gLastPhase[MAX_FRAME_LENGTH/2+1];
    double gSumPhase[MAX_FRAME_LENGTH/2+1];
    double gOutputAccum[2*MAX_FRAME_LENGTH];
    double gAnaFreq[MAX_FRAME_LENGTH];
    double gAnaMagn[MAX_FRAME_LENGTH];
    double gSynFreq[MAX_FRAME_LENGTH];
    double gSynMagn[MAX_FRAME_LENGTH];

    double freqPerBin, expct;
    long inFifoLatency, fftFrameSize2;

    long stepSize;
    long osamp;
    double sampleRate;
    long fftFrameSize;
    double window[MAX_FRAME_LENGTH];
    fftw_plan pfft;
    fftw_plan pifft;
};

smbStruct smbData;


void smbPitchShiftInit(long stepSize_i, long fftFrameSize_i, double sampleRate_i){
    
    
    smbData.stepSize = stepSize_i;
    smbData.osamp = fftFrameSize_i/stepSize_i;
    smbData.sampleRate = sampleRate_i;
    
    smbData.fftFrameSize = fftFrameSize_i;
	smbData.fftFrameSize2 = smbData.fftFrameSize/2;
	smbData.freqPerBin = smbData.sampleRate/(double)smbData.fftFrameSize;
	smbData.expct = 2.*M_PI*(double)smbData.stepSize/(double)smbData.fftFrameSize;
	smbData.inFifoLatency = smbData.fftFrameSize-smbData.stepSize;
    
    memset(smbData.gInFIFO, 0, MAX_FRAME_LENGTH*sizeof(double));
    memset(smbData.gFFTworksp, 0, 2*MAX_FRAME_LENGTH*sizeof(double));
    memset(smbData.gLastPhase, 0, (MAX_FRAME_LENGTH/2+1)*sizeof(double));
    memset(smbData.gSumPhase, 0, (MAX_FRAME_LENGTH/2+1)*sizeof(double));
    memset(smbData.gOutputAccum, 0, 2*MAX_FRAME_LENGTH*sizeof(double));
    memset(smbData.gAnaFreq, 0, MAX_FRAME_LENGTH*sizeof(double));
    memset(smbData.gAnaMagn, 0, MAX_FRAME_LENGTH*sizeof(double));

    for(int k = 0; k < smbData.fftFrameSize;++k){
        smbData.window[k] = -.5*cos(2.*M_PI*(double)k/(double)smbData.fftFrameSize)+.5;
    }

    //fftw_destroy_plan(pfft);
    //fftw_destroy_plan(pfft);
    //pfft = fftw_plan_dft_1d(fftFrameSize, (fftw_complex *)gFFTworksp, (fftwf_complex *)gFFTworksp, FFTW_FORWARD, FFTW_ESTIMATE);
    //pifft = fftw_plan_dft_1d(fftFrameSize, (fftw_complex *)gFFTworksp, (fftwf_complex *)gFFTworksp, FFTW_BACKWARD, FFTW_ESTIMATE);

    smbData.pfft = fftw_plan_dft_r2c_1d(smbData.fftFrameSize, smbData.gFFTworksp_real, (fftw_complex *)smbData.gFFTworksp,FFTW_MEASURE);
    smbData.pifft = fftw_plan_dft_c2r_1d(smbData.fftFrameSize, (fftw_complex *)smbData.gFFTworksp, smbData.gFFTworksp_real,FFTW_MEASURE);
}

void smbPitchShift(double pitchShift, double *indata, double *outdata)
/*
	Routine smbPitchShift(). See top of file for explanation
	Purpose: doing pitch shifting while maintaining duration using the Short
	Time Fourier Transform.
	Author: (c)1999-2015 Stephan M. Bernsee <s.bernsee [AT] zynaptiq [DOT] com>
*/
{

    double magn, phase, tmp, real, imag;
    long k, qpd, index;
    
    memcpy(&smbData.gInFIFO[smbData.inFifoLatency], indata, (size_t) smbData.stepSize*sizeof(double));
    
    /* do windowing and re,im interleave */
    for (k = 0; k < smbData.fftFrameSize;k++) {
        smbData.gFFTworksp_real[k] = smbData.window[k]*smbData.gInFIFO[k];
//         gFFTworksp[2*k] = gInFIFO[k] * window[k];
//         gFFTworksp[2*k+1] = 0.;

    }

    /* ***************** ANALYSIS ******************* */
    /* do transform */
    fftw_execute(smbData.pfft);
    
    /* this is the analysis step */
    for (k = 0; k <= smbData.fftFrameSize2; k++) {

        /* de-interlace FFT buffer */
        real = smbData.gFFTworksp[2*k];
        imag = smbData.gFFTworksp[2*k+1];

        /* compute magnitude and phase */
        magn = sqrt(real*real + imag*imag);
        phase = atan2(imag,real);

        /* compute phase difference */
        tmp = (phase - smbData.gLastPhase[k]);
        smbData.gLastPhase[k] = phase;

        /* subtract expected phase difference */
        tmp -= (double)k*smbData.expct;

        /* map delta phase into +/- Pi interval */
        qpd = tmp/M_PI;
        if (qpd >= 0) qpd += qpd&1;
        else qpd -= qpd&1;
        tmp -= M_PI*(double)qpd;

        /* get deviation from bin frequency from the +/- Pi interval */
        tmp = smbData.osamp*tmp/(2.*M_PI);

        /* compute the k-th partials' true frequency */
        tmp = (double)k*smbData.freqPerBin + tmp*smbData.freqPerBin;

        /* store magnitude and true frequency in analysis arrays */
        smbData.gAnaMagn[k] = magn;
        smbData.gAnaFreq[k] = tmp;

    }

    /* ***************** PROCESSING ******************* */
    /* this does the actual pitch shifting */
    memset(smbData.gSynMagn, 0, smbData.fftFrameSize*sizeof(double));
    memset(smbData.gSynFreq, 0, smbData.fftFrameSize*sizeof(double));
    for (k = 0; k <= smbData.fftFrameSize2; k++) { 
        index = k*pitchShift;
        if (index <= smbData.fftFrameSize2) { 
            smbData.gSynMagn[index] += smbData.gAnaMagn[k]; 
            smbData.gSynFreq[index] = smbData.gAnaFreq[k] * pitchShift; 
        } 
    }

    /* ***************** SYNTHESIS ******************* */
    /* this is the synthesis step */
    for (k = 0; k <= smbData.fftFrameSize2; k++) {

        /* get magnitude and true frequency from synthesis arrays */
        magn = smbData.gSynMagn[k];
        tmp = smbData.gSynFreq[k];

        /* subtract bin mid frequency */
        tmp -= (double)k*smbData.freqPerBin;

        /* get bin deviation from freq deviation */
        tmp /= smbData.freqPerBin;

        /* take osamp into account */
        tmp = 2.*M_PI*tmp/smbData.osamp;

        /* add the overlap phase advance back in */
        tmp += (double)k*smbData.expct;

        /* accumulate delta phase to get bin phase */
        smbData.gSumPhase[k] += tmp;
        phase = smbData.gSumPhase[k];

        /* get real and imag part and re-interleave */
        smbData.gFFTworksp[2*k] = magn*cos(phase);
        smbData.gFFTworksp[2*k+1] = magn*sin(phase);
    } 

    /* zero negative frequencies */
    //for (k = fftFrameSize+2; k < 2*fftFrameSize; k++) gFFTworksp[k] = 0.;

    /* do inverse transform */
    fftw_execute(smbData.pifft);
    
    for(k=0; k < smbData.fftFrameSize; k++) {
        smbData.gFFTworksp_real[k] = smbData.window[k]*smbData.gFFTworksp_real[k]/smbData.fftFrameSize/(smbData.osamp*3/(double)8);
    }

    /* do windowing and add to output accumulator */ 
    for(k=0; k < smbData.fftFrameSize; k++) {
        smbData.gOutputAccum[k] += smbData.gFFTworksp_real[k];
    }
    for (k = 0; k < smbData.stepSize; k++) outdata[k] = smbData.gOutputAccum[k];

    /* shift accumulator */
    memmove(smbData.gOutputAccum, smbData.gOutputAccum+smbData.stepSize, smbData.fftFrameSize*sizeof(double));

    /* move input FIFO */
    for (k = 0; k < smbData.inFifoLatency; k++) smbData.gInFIFO[k] = smbData.gInFIFO[k+smbData.stepSize];
    
    return;
}

