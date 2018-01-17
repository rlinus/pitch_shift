#include <iostream>
#include <stdlib.h>
#include <string>
#include <math.h>
#include "mex.h"
#include "Stk.h"
#include "RtAudio.h"

#include "cpvPitchShift.h"
#include "smbPitchShift.h"

#define RUBBERBAND //if RUBBERBAND lib is available

#if defined RUBBERBAND
#include "RubberBandStretcher.h"
using namespace RubberBand;
#endif


using namespace std;
using namespace stk;

void volume_normalizer(double *frame);

bool volume_normalization = false;
double volume_normalization_timeconst = 0;
double volume_normalization_delta = 0;

int deviceId = 0;
int shifterId = 0;
int windowSize = 1024;

const int frameSize = 64;
const int sampleRate = 44100;

RtAudio dac;



#if defined RUBBERBAND
    RubberBandStretcher::Options options =  RubberBandStretcher::OptionProcessRealTime |
                                        RubberBandStretcher::OptionFormantPreserved |
                                        RubberBandStretcher::OptionWindowShort |
                                        //RubberBandStretcher::OptionSmoothingOn |
                                        //RubberBandStretcher::OptionPhaseIndependent |
                                        RubberBandStretcher::OptionTransientsSmooth |
                                        RubberBandStretcher::OptionPitchHighQuality;//Consistency;
    RubberBandStretcher shifter(sampleRate,1,options,1.0,1.0);
#endif


double output[frameSize];
double *output_p = output;
int i_frame;

const size_t data_array_length = 20*sampleRate;
double signal[data_array_length];



double pitch_factor;
double prev_pitch_factor = 0;
int signal_length;


void init(void){
    Stk::setSampleRate(sampleRate);
    
    i_frame = -1;
    
    pitch_factor = 1;
    
    switch(shifterId){
        case 0:
            smbPitchShiftInit(frameSize, windowSize, sampleRate);
            break;
        case 1:
            cpvPitchShiftInit(frameSize, windowSize, sampleRate);
            break;
        #if defined RUBBERBAND
        case 2:
            shifter.reset();
            break;
        #endif   
    }
}

int tick( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
         double streamTime, RtAudioStreamStatus status, void *dataPointer )
{
    i_frame++;

// 	if(status == RTAUDIO_INPUT_OVERFLOW){
// 		mexPrintf("overflow!\n");
// 	} else if (status == RTAUDIO_OUTPUT_UNDERFLOW){
// 		mexPrintf("underflow!\n");
// 	}
    
    if(status) mexPrintf("underflow/overflow detected!\n");
    
    if(nBufferFrames!=frameSize){
        mexPrintf("not enough samples received: %i,frame: %i\n",nBufferFrames, i_frame);
    }

	double *ibuffer = &signal[i_frame*frameSize%signal_length];
	double *obuffer = (double *) outputBuffer;
    
    
    switch(shifterId){
        case 0:
            smbPitchShift(pitch_factor, ibuffer, output);
            break;
        case 1:
            cpvPitchShift(pitch_factor, ibuffer, output);
            break;
        #if defined RUBBERBAND
        case 2:
            float in[frameSize];
            float out[frameSize];
            float * in_p = in;
            float * out_p = out;
            for(int i=0; i < frameSize; ++i){
                in[i] = (float) ibuffer[i];
            }

            shifter.setPitchScale(pitch_factor);
            shifter.process(&in_p, frameSize,0);
            int r = shifter.available();

            if(r>=frameSize && i_frame+1>=18){
                shifter.retrieve(&out_p, frameSize);
                for(int i=0; i < frameSize; ++i){
                    output[i] = out[i];
                }
            }else{
                for(int i=0; i < frameSize; ++i){
                    output[i] = 0;
                }
            } 
            break;
        #endif
    }

	//normalize volume
	if (volume_normalization) volume_normalizer(output);
    
    if(prev_pitch_factor == pitch_factor){    
        for (int i=0; i<frameSize; i++ ){  
                double out = output[i];

                if(out>1.0f) out = 1.0f;
                if(out<-1.0f) out = -1.0f;
                *obuffer++ = out;
                *obuffer++ = out;
        }
    } else {
        for (int i=0; i<frameSize; i++ ){  
                *obuffer++ = 0;
                *obuffer++ = 0;
        }
    }
    prev_pitch_factor = pitch_factor;
	return 0;
}

void start_stream(void)
{   
    
    RtAudio::StreamParameters oparameters;
    oparameters.deviceId = deviceId; //dac.getDefaultOutputDevice();
	oparameters.nChannels = 2;

	RtAudio::StreamParameters iparameters;
	iparameters.deviceId = deviceId; //dac.getDefaultInputDevice();
	iparameters.nChannels = 0;
    
    RtAudio::StreamOptions options;
    options.flags = RTAUDIO_HOG_DEVICE;// | RTAUDIO_SCHEDULE_REALTIME | RTAUDIO_MINIMIZE_LATENCY;
    //options.priority = 10000000;
    
    unsigned int bufferFrames = frameSize;
    
	try {
		dac.openStream( &oparameters, NULL, RTAUDIO_FLOAT64, (unsigned int)Stk::sampleRate(), &bufferFrames, &tick, NULL, &options);
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
    
    //mexPrintf("HW latency: %i\n",lat_hw);
    
    
    return;
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
//     mexPrintf("N: %i\n",dac.getDeviceCount());
//     RtAudio::DeviceInfo dI = dac.getDeviceInfo(0);
//     string a = "hallo";//dI.name;
//     mexPrintf("P: %i, N: %s, oC: %i, iC: %i\n",dI.probed,a,dI.outputChannels,dI.inputChannels);
//     return;
    if(nrhs == 0){
        if(nlhs>0) {
            mexErrMsgIdAndTxt("rt_pitch_shifter:nlhs", "Too many output arguments.");
        }
        return;
    } else {
        if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfElements(prhs[0])!=1) {
            mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "First input must be a scalar.");
            return;
        }
        
        int v1 = (int) floor(mxGetScalar(prhs[0]));
        
        if(v1 < 0){
            try {
                if(dac.isStreamRunning()) dac.stopStream();
                dac.closeStream();
            }
            catch ( RtAudioError &error ) {
                error.printMessage();
                return;
            }

            if(nlhs>0) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:nlhs", "Too many output arguments.");
            }
            
            //mexPrintf("Number of processed frames: %i\n",i_frame+1);

            return;
        } else if(v1==0){
            if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetNumberOfElements(prhs[1])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Second input must be a scalar.");
                return;
            }
            pitch_factor = mxGetScalar(prhs[1]);
            return;
        } else if(v1==1){
            
            if(dac.isStreamRunning() || dac.isStreamOpen()){
                mexErrMsgIdAndTxt("rt_pitch_shifter:openStream", "Stream is already open.");
                return;
            }   
            
            if( nrhs < 2 || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetM(prhs[1])>1 && mxGetN(prhs[1])>1)) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notVector", "Second input must be a vector.");
                return;
            }
            
            signal_length = ((int)(mxGetNumberOfElements(prhs[1])/frameSize))*frameSize;

            double *prhs1 = mxGetPr(prhs[1]);
            
            for(int i = 0; i < signal_length; ++i){
                signal[i] = prhs1[i];
            }
            
            pitch_factor = 1;
            shifterId = 0;
            deviceId = 0;
			volume_normalization = false;
			volume_normalization_timeconst = 0.5;
			volume_normalization_delta = 0.5;
            
            if(nrhs>2) pitch_factor = mxGetScalar(prhs[2]);
            if(nrhs>3) shifterId =  mxGetScalar(prhs[3]);
            if(nrhs>4) deviceId = mxGetScalar(prhs[4]);
			if(nrhs>5) volume_normalization = (bool) mxGetScalar(prhs[5]);
			if (nrhs>6) volume_normalization_timeconst = mxGetScalar(prhs[6]);
			if (volume_normalization_timeconst < 0) volume_normalization_timeconst = 0;
			if (volume_normalization_timeconst > 1) volume_normalization_timeconst = 1;

			if (nrhs>7) volume_normalization_delta = mxGetScalar(prhs[7]);
			if (volume_normalization_delta <= 0) volume_normalization_delta = 1e-10;

			if (nrhs>8) {
				int tid = int(mxGetScalar(prhs[8]));
				if (tid == 0) {
					windowSize = 512;
				}
				else if (tid == 2) {
					windowSize = 2048;
				}
				else {
					windowSize = 1024;
				}
			}
			else {
				windowSize = 1024;
			}
            
            init();
            start_stream();
            
        } else {
            mexErrMsgIdAndTxt("rt_pitch_shifter:invalidInput", "Invalid input");
        }

    }
    
    return;
}

//smooth power of output signal
void volume_normalizer(double *frame) {
	static double filterd_power = 0;

	for (int i = 0; i < frameSize; ++i) {
		filterd_power = (1 - volume_normalization_timeconst)*filterd_power + volume_normalization_timeconst * frame[i] * frame[i];
		double x = sqrt(filterd_power) / volume_normalization_delta;
		if (x>0) frame[i] = frame[i] * log(x + 1) / x;
	}
}

