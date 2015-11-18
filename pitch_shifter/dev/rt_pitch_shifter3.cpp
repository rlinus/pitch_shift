#include <iostream>
#include <stdlib.h>
#include <random>
#include "mex.h"
#include "Stk.h"
#include "RtAudio.h"
#include "RubberBandStretcher.h"

using namespace std;
using namespace stk;
using namespace RubberBand;

float frameAmplitude(float *samples);

const int frameSize = 64;
const int sampleRate = 44100;                       

float *output;
int i_frame;

int i_low_amp;
const float threshold = 0.005f;

const size_t data_array_length = 360*44100;
float input_signal[data_array_length];
float output_signal[data_array_length];

RtAudio dac;

RubberBandStretcher::Options options =  RubberBandStretcher::OptionProcessRealTime |
                                        RubberBandStretcher::OptionFormantPreserved |
                                        RubberBandStretcher::OptionWindowShort |
                                        RubberBandStretcher::OptionPitchHighConsistency;
RubberBandStretcher shifter(sampleRate,1,options,1.0,pow(2.0,0/12.0));

double pitch_factors[1000] ;
int pitchchange_times_f[1000];
int num_changes = 0;
int i_change = 0;



void init_global_vars(void){
    i_frame = -1;
    i_low_amp = 0;
    
    i_change = 0;
    
    shifter.reset();
    shifter.setPitchScale(1.0);
    
    memset(input_signal,0,sizeof(input_signal));
    memset(output_signal,0,sizeof(output_signal));
}

int tick( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
         double streamTime, RtAudioStreamStatus status, void *dataPointer )
{
    i_frame++;

	if(status == RTAUDIO_INPUT_OVERFLOW){
		mexPrintf("overflow!\n");
	} else if (status == RTAUDIO_OUTPUT_UNDERFLOW){
		mexPrintf("underflow!\n");
	}
    
    if(nBufferFrames!=frameSize){
        mexPrintf("not enough samples received: %i,frame: %i\n",nBufferFrames, i_frame);
    }

	register float *ibuffer = (float *) inputBuffer;
	register float *obuffer = (float *) outputBuffer;
    
    memcpy(&input_signal[i_frame*frameSize], ibuffer,sizeof(input_signal[0])*nBufferFrames);
    
    float amp = frameAmplitude(ibuffer);
    
    if(amp < threshold){
        ++i_low_amp;
    }else{
        i_low_amp = 0;
    }
    
    if(i_frame == pitchchange_times_f[i_change]){
        shifter.setPitchScale(pitch_factors[i_change]);
        if(i_change < num_changes-1) i_change++;
    }
     
    shifter.process(&ibuffer, nBufferFrames,0);
    
    int r = shifter.retrieve(&output, nBufferFrames);
    
    memcpy(&output_signal[i_frame*frameSize],output,r*sizeof(float));
    

	for (int i=0; i<frameSize; i++ ){       
		*obuffer++ = output_signal[i_frame*frameSize+i];
    	*obuffer++ = output_signal[i_frame*frameSize+i];
	}
    
	return 0;
}


void rt_pitch_shifter(void)
{    
    
    output = new float[frameSize];
    
    Stk::setSampleRate( sampleRate );
    
    RtAudio::StreamParameters oparameters;
    oparameters.deviceId = 1; //dac.getDefaultOutputDevice();
	oparameters.nChannels = 2;

	RtAudio::StreamParameters iparameters;
	iparameters.deviceId = 1; //dac.getDefaultInputDevice();
	iparameters.nChannels = 1;
    
    unsigned int bufferFrames = frameSize;
    
    
    shifter.setMaxProcessSize(bufferFrames);
    
    
	try {
		dac.openStream( &oparameters, &iparameters, RTAUDIO_FLOAT32, (unsigned int)Stk::sampleRate(), &bufferFrames, &tick, NULL, NULL);
	}
	catch ( RtAudioError &error ) {
		error.printMessage();
		return;
	}
    
    try {
		dac.startStream();
	}
	catch ( RtAudioError &error ) {
		error.printMessage();
		return;
	}
    
    int lat_hw = dac.getStreamLatency();
    int lat_sw = shifter.getLatency();
    
    mexPrintf("HW latency: %i, SW latency: %i\n",lat_hw,lat_sw);
    
    
    return;
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
    
    
    if(nrhs == 0){
        if(nlhs>0) {
            mexErrMsgIdAndTxt("rt_pitch_shifter:nlhs", "Too many output arguments.");
        }
        return;
    } else {
        if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || (mxGetM(prhs[0])>1 && mxGetN(prhs[0])>1)) {
            mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notVector", "First input must be a vector.");
            return;
        }
        
        int v1 = (int) floor(mxGetScalar(prhs[0]));
        
        if(v1 < 0){
            try {
                dac.closeStream();
            }
            catch ( RtAudioError &error ) {
                error.printMessage();
                return;
            }

            delete[] output;

            mxArray *input_signal_m = mxCreateDoubleMatrix((i_frame+1)*frameSize,1,mxREAL);
            double *input_signal_mp = mxGetPr(input_signal_m);
            mxArray *output_signal_m = mxCreateDoubleMatrix((i_frame+1)*frameSize,1,mxREAL);
            double *output_signal_mp = mxGetPr(output_signal_m);

            for(int i = 0; i<(i_frame+1)*frameSize; ++i){
                input_signal_mp[i] = (double) input_signal[i];
                output_signal_mp[i] = (double) output_signal[i];
            }


            if(nlhs >= 1) plhs[0] = input_signal_m;
            if(nlhs >= 2) plhs[1] = output_signal_m;
            
            mexPrintf("Number of processed frames: %i\n",i_frame);

            return;
        } else {
            if(nlhs>0) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:nlhs", "Too many output arguments.");
            }
            
            if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetM(prhs[0])>1 && mxGetN(prhs[0])>1)) {
                mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notVector", "Second input must be a vector.");
                return;
            }
            
            if(mxGetNumberOfElements(prhs[0]) != mxGetNumberOfElements(prhs[1])){
                mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notSameLengths", "The input vectors must have the same length.");
                return;
            }
            
            num_changes = mxGetNumberOfElements(prhs[0]);
            
            double *prhs0 = mxGetPr(prhs[0]);
            double *prhs1 = mxGetPr(prhs[1]);
            
            for(int i = 0; i < num_changes; ++i){
                pitchchange_times_f[i] = (int) floor(prhs0[i]);
                pitch_factors[i] = prhs1[i];
            }
            
            init_global_vars();
            rt_pitch_shifter();
        }
    }
    
    return;
}

float frameAmplitude(float *samples) {
    float sum = 0;
    for(int i = 0; i < frameSize; ++i) {
        sum += abs(samples[0]);
    }
    return sum/frameSize;
        
}