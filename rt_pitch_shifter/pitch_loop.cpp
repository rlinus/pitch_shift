#include <iostream>
#include <stdlib.h>
#include <random>
#include <string>
#include <math.h>
#include "mex.h"
#include "Stk.h"
#include "RtAudio.h"
#include "Wurley.h"
#include "Rhodey.h"
#include "BeeThree.h"
#include "Drummer.h"
#include "dywapitchtrack.h"



#if defined CPVSHIFTER
    #include "cpvPitchShift.h"
#elif defined SMBSHIFTER
    #include "smbPitchShift.h"
#elif defined DIRACSHIFTER
    #include "Dirac.h"
#elif defined RUBBERBAND
    #include "RubberBandStretcher.h"
    using namespace RubberBand;
#else
    #error "PITCHSHIFTER is not defined!"
#endif

#define DEVICEID 0
//#define RAWWAVEPATH "C:/Users/Linus/Documents/MATLAB/pitch_shift/rt_pitch_shifter/stk-4.5.0/rawwaves/"


using namespace std;
using namespace stk;



const int frameSize = 64;
const int sampleRate = 44100;

RtAudio dac;


#if defined DIRACSHIFTER
    void *diracFx = DiracFxCreate(kDiracQualityBest, sampleRate, 1);
#elif defined RUBBERBAND
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

const size_t data_array_length = 10*sampleRate;
double signal[data_array_length];



double pitch_factor;
int signal_length;


void init(void){
    Stk::setSampleRate(sampleRate);
    Stk::setRawwavePath(RAWWAVEPATH);
    
    i_frame = -1;
    
    pitch_factor = 1;
    
    #if defined CPVSHIFTER
    	cpvPitchShiftInit(frameSize, 1024/frameSize, sampleRate);
    #elif defined SMBSHIFTER
        smbPitchShiftInit(frameSize, 1024/frameSize, sampleRate);
    #elif defined RUBBERBAND
        shifter.reset();
    #endif
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
    
    
    #if defined CPVSHIFTER
    	cpvPitchShift(pitch_factor, ibuffer, output);
    #elif defined SMBSHIFTER
        smbPitchShift(pitch_factor, output);
    #elif defined DIRACSHIFTER
        float in[frameSize];
        float out[frameSize];
        for(int i=0; i < frameSize; ++i){
            in[i] = (float) ibuffer[i];
        }
        
        DiracFxProcessFloatInterleaved(1.0, pitch_factor, in, out, frameSize, diracFx);
        for(int i=0; i < frameSize; ++i){
            output[i] = out[i];
        }
    #elif defined RUBBERBAND
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
    #endif
    
        
    for (int i=0; i<frameSize; i++ ){  
            double out = output[i];
            
            if(out>1.0f) out = 1.0f;
            if(out<-1.0f) out = -1.0f;
            *obuffer++ = out;
            *obuffer++ = out;
    }
	return 0;
}

void start_stream(void)
{   
    
    RtAudio::StreamParameters oparameters;
    oparameters.deviceId = DEVICEID; //dac.getDefaultOutputDevice();
	oparameters.nChannels = 2;

	RtAudio::StreamParameters iparameters;
	iparameters.deviceId = DEVICEID; //dac.getDefaultInputDevice();
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
    
    #ifdef DIRACSHIFTER
        mexPrintf("SW latency: %i\n",DiracFxLatencyFrames(sampleRate));
    #endif
    
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
            
            if( nrhs < 2 || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetM(prhs[0])>1 && mxGetN(prhs[0])>1)) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notVector", "Second input must be a vector.");
                return;
            }
            
            signal_length = ((int)(mxGetNumberOfElements(prhs[1])/frameSize))*frameSize;

            double *prhs1 = mxGetPr(prhs[1]);
            
            for(int i = 0; i < signal_length; ++i){
                signal[i] = prhs1[i];
            }
            
            init();
            start_stream();
            
        } else {
            mexErrMsgIdAndTxt("rt_pitch_shifter:invalidInput", "Invalid input");
        }

    }
    
    return;
}

