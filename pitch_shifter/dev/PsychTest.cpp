#include <iostream>
#include <stdlib.h>
#include <string>
#include <math.h>
#include "mex.h"
#include "Stk.h"
#include "RtAudio.h"




using namespace std;
using namespace stk;

int deviceId = 0;


std::string get_mex_path(void);


const int frameSize = 64;
const int sampleRate = 44100;

const size_t data_array_length = 120*sampleRate;
double input_signal[data_array_length];
double output_signal[data_array_length];


RtAudio dac;

double output[frameSize];
double *output_p = output;
int i_frame;

int is_finished = 0;





void init(void){
    Stk::setSampleRate(sampleRate);
    
    i_frame = -1;
    
    is_finished = 0;

    memset(input_signal,0,sizeof(input_signal));
    memset(output_signal,0,sizeof(output_signal));


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
    //mexPrintf("st: %i\n", status);
    
    if(status) mexPrintf("underflow/overflow detected!\n");
    
    if(nBufferFrames!=frameSize){
        mexPrintf("not enough samples received: %i,frame: %i\n",nBufferFrames, i_frame);
    }

	double *ibuffer = (double *) inputBuffer;
	double *obuffer = (double *) outputBuffer;
    
    
    //ibuffer = &beep_sound[i_frame*frameSize];
    
    memcpy(&input_signal[i_frame*frameSize], ibuffer,sizeof(input_signal[0])*nBufferFrames);
    

        
    is_finished = 1;
    
        
    for (int i=0; i<nBufferFrames; i++ ){  

            *obuffer++ = ibuffer[i];
            *obuffer++ = ibuffer[i];
//             double tem = sin(1000*(i_frame+i)/sampleRate);
//             *obuffer++ = tem;
//             *obuffer++ = tem;
    }
	return 0;
}

void errorcb(RtAudioError::Type type, const std::string &errorText){
    mexPrintf("ERROR: %i\n",type);
}


void start_stream(void)
{   
    
    RtAudio::StreamParameters oparameters;
    oparameters.deviceId = deviceId; //dac.getDefaultOutputDevice();
	oparameters.nChannels = 2;

	RtAudio::StreamParameters iparameters;
	iparameters.deviceId = deviceId; //dac.getDefaultInputDevice();
	iparameters.nChannels = 1;
    
    RtAudio::StreamOptions options;
    options.flags = RTAUDIO_HOG_DEVICE;// | RTAUDIO_SCHEDULE_REALTIME;
    //options.priority = 10000000;
    
    unsigned int bufferFrames = frameSize;
    
    
	try {
		dac.openStream( &oparameters, &iparameters, RTAUDIO_FLOAT64, sampleRate, &bufferFrames, &tick, NULL, &options, &errorcb);
        mexPrintf("actual frameSize: %i\n",bufferFrames);
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
    
    mexPrintf("HW latency: %i\n",lat_hw);

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
            
            //mexPrintf("Number of processed frames: %i\n",i_frame+1);

            return;
        } else if(v1==0){
            plhs[0] = mxCreateDoubleScalar(is_finished);
            return;
        } else {
            
            if(dac.isStreamRunning() || dac.isStreamOpen()){
                mexErrMsgIdAndTxt("rt_pitch_shifter:openStream", "Stream is already open.");
                return;
            }     
            
            if(nlhs>0) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:nlhs", "Too many output arguments.");
                return;
            }

     
            
            
            //start stream
            init();
            start_stream();
        }
    }
    
    return;
}


std::string get_mex_path(void){
    mxArray *rhs[1], *lhs[1];
    char *path, *name;
    size_t lenpath, lenname, n;
    rhs[0] = mxCreateString("fullpath");
    mexCallMATLAB(1, lhs, 1, rhs, "mfilename");
    mxDestroyArray(rhs[0]);
    path = mxArrayToString(lhs[0]);
    mxDestroyArray(lhs[0]);
    mexCallMATLAB(1, lhs, 0, rhs, "mfilename");
    name = mxArrayToString(lhs[0]);
    mxDestroyArray(lhs[0]);

    lenpath = strlen(path);
    lenname = strlen(name);
    n = lenpath - lenname;
    path[n] = '\0';
    
    
    for(int i=0; i<n; ++i){
        if(path[i]=='\\') path[i]='/';
    }
    
    
    std::string str(path);
    return str;
}