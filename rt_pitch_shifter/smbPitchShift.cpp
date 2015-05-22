/****************************************************************************
*
* NAME: smbPitchShift.cpp
* VERSION: 1.2
* HOME URL: http://blogs.zynaptiq.com/bernsee
* KNOWN BUGS: none
*
* SYNOPSIS: Routine for doing pitch shifting while maintaining
* duration using the Short Time Fourier Transform.
*
* DESCRIPTION: The routine takes a pitchShift factor value which is between 0.5
* (one octave down) and 2. (one octave up). A value of exactly 1 does not change
* the pitch. numSampsToProcess tells the routine how many samples in indata[0...
* numSampsToProcess-1] should be pitch shifted and moved to outdata[0 ...
* numSampsToProcess-1]. The two buffers can be identical (ie. it can process the
* data in-place). fftFrameSize defines the FFT frame size used for the
* processing. Typical values are 1024, 2048 and 4096. It may be any value <=
* MAX_FRAME_LENGTH but it MUST be a power of 2. osamp is the STFT
* oversampling factor which also determines the overlap between adjacent STFT
* frames. It should at least be 4 for moderate scaling ratios. A value of 32 is
* recommended for best quality. sampleRate takes the sample rate for the signal 
* in unit Hz, ie. 44100 for 44.1 kHz audio. The data passed to the routine in 
* indata[] should be in the range [-1.0, 1.0), which is also the output range 
* for the data, make sure you scale the data accordingly (for 16bit signed integers
* you would have to divide (and multiply) by 32768). 
*
* COPYRIGHT 1999-2015 Stephan M. Bernsee <s.bernsee [AT] zynaptiq [DOT] com>
*
* 						The Wide Open License (WOL)
*
* Permission to use, copy, modify, distribute and sell this software and its
* documentation for any purpose is hereby granted without fee, provided that
* the above copyright notice and this license appear in all source copies. 
* THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY OF
* ANY KIND. See http://www.dspguru.com/wol.htm for more information.
*
*****************************************************************************/ 
#include "smbPitchShift.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "fftw3.h"


#define M_PI 3.14159265358979323846
#define MAX_FRAME_LENGTH 8192


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


// -----------------------------------------------------------------------------------------------------------------


void smbPitchShiftInit(long stepSize_i, long osamp_i, double sampleRate_i){
    stepSize = stepSize_i;
    osamp = osamp_i;
    sampleRate = sampleRate_i;
    
    fftFrameSize = stepSize*osamp;
	fftFrameSize2 = fftFrameSize/2;
	freqPerBin = sampleRate/(double)fftFrameSize;
	expct = 2.*M_PI*(double)stepSize/(double)fftFrameSize;
	inFifoLatency = fftFrameSize-stepSize;
    
    memset(gInFIFO, 0, MAX_FRAME_LENGTH*sizeof(double));
    memset(gFFTworksp, 0, 2*MAX_FRAME_LENGTH*sizeof(double));
    memset(gLastPhase, 0, (MAX_FRAME_LENGTH/2+1)*sizeof(double));
    memset(gSumPhase, 0, (MAX_FRAME_LENGTH/2+1)*sizeof(double));
    memset(gOutputAccum, 0, 2*MAX_FRAME_LENGTH*sizeof(double));
    memset(gAnaFreq, 0, MAX_FRAME_LENGTH*sizeof(double));
    memset(gAnaMagn, 0, MAX_FRAME_LENGTH*sizeof(double));

    for(int k = 0; k < fftFrameSize;++k){
        window[k] = -.5*cos(2.*M_PI*(double)k/(double)fftFrameSize)+.5;
    }

    //fftw_destroy_plan(pfft);
    //fftw_destroy_plan(pfft);
    //pfft = fftw_plan_dft_1d(fftFrameSize, (fftw_complex *)gFFTworksp, (fftwf_complex *)gFFTworksp, FFTW_FORWARD, FFTW_ESTIMATE);
    //pifft = fftw_plan_dft_1d(fftFrameSize, (fftw_complex *)gFFTworksp, (fftwf_complex *)gFFTworksp, FFTW_BACKWARD, FFTW_ESTIMATE);

    pfft = fftw_plan_dft_r2c_1d(fftFrameSize, gFFTworksp_real, (fftw_complex *)gFFTworksp,FFTW_MEASURE);
    pifft = fftw_plan_dft_c2r_1d(fftFrameSize, (fftw_complex *)gFFTworksp, gFFTworksp_real,FFTW_MEASURE);
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
    
    memcpy(&gInFIFO[inFifoLatency], indata, (size_t) stepSize*sizeof(double));
    
    /* do windowing and re,im interleave */
    for (k = 0; k < fftFrameSize;k++) {
        gFFTworksp_real[k] = window[k]*gInFIFO[k];
//         gFFTworksp[2*k] = gInFIFO[k] * window[k];
//         gFFTworksp[2*k+1] = 0.;

    }

    /* ***************** ANALYSIS ******************* */
    /* do transform */
    fftw_execute(pfft);
    
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

        /* map delta phase into +/- Pi interval */
        qpd = tmp/M_PI;
        if (qpd >= 0) qpd += qpd&1;
        else qpd -= qpd&1;
        tmp -= M_PI*(double)qpd;

        /* get deviation from bin frequency from the +/- Pi interval */
        tmp = osamp*tmp/(2.*M_PI);

        /* compute the k-th partials' true frequency */
        tmp = (double)k*freqPerBin + tmp*freqPerBin;

        /* store magnitude and true frequency in analysis arrays */
        gAnaMagn[k] = magn;
        gAnaFreq[k] = tmp;

    }

    /* ***************** PROCESSING ******************* */
    /* this does the actual pitch shifting */
    memset(gSynMagn, 0, fftFrameSize*sizeof(double));
    memset(gSynFreq, 0, fftFrameSize*sizeof(double));
    for (k = 0; k <= fftFrameSize2; k++) { 
        index = k*pitchShift;
        if (index <= fftFrameSize2) { 
            gSynMagn[index] += gAnaMagn[k]; 
            gSynFreq[index] = gAnaFreq[k] * pitchShift; 
        } 
    }

    /* ***************** SYNTHESIS ******************* */
    /* this is the synthesis step */
    for (k = 0; k <= fftFrameSize2; k++) {

        /* get magnitude and true frequency from synthesis arrays */
        magn = gSynMagn[k];
        tmp = gSynFreq[k];

        /* subtract bin mid frequency */
        tmp -= (double)k*freqPerBin;

        /* get bin deviation from freq deviation */
        tmp /= freqPerBin;

        /* take osamp into account */
        tmp = 2.*M_PI*tmp/osamp;

        /* add the overlap phase advance back in */
        tmp += (double)k*expct;

        /* accumulate delta phase to get bin phase */
        gSumPhase[k] += tmp;
        phase = gSumPhase[k];

        /* get real and imag part and re-interleave */
        gFFTworksp[2*k] = magn*cos(phase);
        gFFTworksp[2*k+1] = magn*sin(phase);
    } 

    /* zero negative frequencies */
    //for (k = fftFrameSize+2; k < 2*fftFrameSize; k++) gFFTworksp[k] = 0.;

    /* do inverse transform */
    fftw_execute(pifft);
    
    for(k=0; k < fftFrameSize; k++) {
        gFFTworksp_real[k] = window[k]*gFFTworksp_real[k]/fftFrameSize/(osamp*3/(double)8);
    }

    /* do windowing and add to output accumulator */ 
    for(k=0; k < fftFrameSize; k++) {
        gOutputAccum[k] += gFFTworksp_real[k];
    }
    for (k = 0; k < stepSize; k++) outdata[k] = gOutputAccum[k];

    /* shift accumulator */
    memmove(gOutputAccum, gOutputAccum+stepSize, fftFrameSize*sizeof(double));

    /* move input FIFO */
    for (k = 0; k < inFifoLatency; k++) gInFIFO[k] = gInFIFO[k+stepSize];
    
    return;
}

// -----------------------------------------------------------------------------------------------------------------


void smbFft(double *fftBuffer, long fftFrameSize, long sign)
/* 
	FFT routine, (C)1996 S.M.Bernsee. Sign = -1 is FFT, 1 is iFFT (inverse)
	Fills fftBuffer[0...2*fftFrameSize-1] with the Fourier transform of the
	time domain data in fftBuffer[0...2*fftFrameSize-1]. The FFT array takes
	and returns the cosine and sine parts in an interleaved manner, ie.
	fftBuffer[0] = cosPart[0], fftBuffer[1] = sinPart[0], asf. fftFrameSize
	must be a power of 2. It expects a complex input signal (see footnote 2),
	ie. when working with 'common' audio signals our input signal has to be
	passed as {in[0],0.,in[1],0.,in[2],0.,...} asf. In that case, the transform
	of the frequencies of interest is in fftBuffer[0...fftFrameSize].
*/
{
	double wr, wi, arg, *p1, *p2, temp;
	double tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
	long i, bitm, j, le, le2, k;

	for (i = 2; i < 2*fftFrameSize-2; i += 2) {
		for (bitm = 2, j = 0; bitm < 2*fftFrameSize; bitm <<= 1) {
			if (i & bitm) j++;
			j <<= 1;
		}
		if (i < j) {
			p1 = fftBuffer+i; p2 = fftBuffer+j;
			temp = *p1; *(p1++) = *p2;
			*(p2++) = temp; temp = *p1;
			*p1 = *p2; *p2 = temp;
		}
	}
	for (k = 0, le = 2; k < (long)(log((double)fftFrameSize)/log((double)2.)+.5); k++) {
		le <<= 1;
		le2 = le>>1;
		ur = 1.0;
		ui = 0.0;
		arg = M_PI / (le2>>1);
		wr = cos(arg);
		wi = sign*sin(arg);
		for (j = 0; j < le2; j += 2) {
			p1r = fftBuffer+j; p1i = p1r+1;
			p2r = p1r+le2; p2i = p2r+1;
			for (i = j; i < 2*fftFrameSize; i += le) {
				tr = *p2r * ur - *p2i * ui;
				ti = *p2r * ui + *p2i * ur;
				*p2r = *p1r - tr; *p2i = *p1i - ti;
				*p1r += tr; *p1i += ti;
				p1r += le; p1i += le;
				p2r += le; p2i += le;
			}
			tr = ur*wr - ui*wi;
			ui = ur*wi + ui*wr;
			ur = tr;
		}
	}
}


// -----------------------------------------------------------------------------------------------------------------

/*

    12/12/02, smb
    
    PLEASE NOTE:
    
    There have been some reports on domain errors when the atan2() function was used
    as in the above code. Usually, a domain error should not interrupt the program flow
    (maybe except in Debug mode) but rather be handled "silently" and a global variable
    should be set according to this error. However, on some occasions people ran into
    this kind of scenario, so a replacement atan2() function is provided here.
    
    If you are experiencing domain errors and your program stops, simply replace all
    instances of atan2() with calls to the smbAtan2() function below.
    
*/


double smbAtan2(double x, double y)
{
  double signx;
  if (x > 0.) signx = 1.;  
  else signx = -1.;
  
  if (x == 0.) return 0.;
  if (y == 0.) return signx * M_PI / 2.;
  
  return atan2(x, y);
}


// -----------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------
