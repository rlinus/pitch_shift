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

float output[frameSize];
float *output_p = output;
int i_frame;
int i_frame_session;

int i_low_amp;
float threshold = 0.005f; //0.05 with intern mic
int minlowampduration_f = 7;

const size_t data_array_length = 10*60*sampleRate;
float input_signal[data_array_length];
float output_signal[data_array_length];

int session_starts[1000];

RtAudio dac;

RubberBandStretcher::Options options =  RubberBandStretcher::OptionProcessRealTime |
                                        RubberBandStretcher::OptionFormantPreserved |
                                        RubberBandStretcher::OptionWindowShort |
                                        RubberBandStretcher::OptionTransientsSmooth |
                                        RubberBandStretcher::OptionPitchHighConsistency;
RubberBandStretcher shifter(sampleRate,1,options,1.0,1.0);

int session_duration_f = 700;
double pitch_factors[1000] ;
int pitchchange_times_f[1000];
int num_sessions = 0;
int i_session = 0;

volatile int is_finished = 0;

bool mark_session_starts = false;

bool add_pink_noise = false;
float noise_gain = 0.1f;
const int num_noise_frames = 5000;
float pink_noise[num_noise_frames*frameSize];
int i_noise_frame = 0;


void init_global_vars(void){
    i_frame_session = session_duration_f;
    i_frame = -1;
    i_low_amp = 0;
    
    is_finished = 0;
    i_session = 0;
    
    i_noise_frame = 0;
    
    shifter.reset();
    shifter.setPitchScale(1.0);
    
    memset(input_signal,0,sizeof(input_signal));
    memset(output_signal,0,sizeof(output_signal));
}

int tick( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
         double streamTime, RtAudioStreamStatus status, void *dataPointer )
{
    i_frame++;
    i_frame_session++;

	if(status == RTAUDIO_INPUT_OVERFLOW){
		mexPrintf("overflow!\n");
	} else if (status == RTAUDIO_OUTPUT_UNDERFLOW){
		mexPrintf("underflow!\n");
	}
    
    if(nBufferFrames!=frameSize){
        mexPrintf("not enough samples received: %i,frame: %i\n",nBufferFrames, i_frame);
    }

	float *ibuffer = (float *) inputBuffer;
	float *obuffer = (float *) outputBuffer;
    
    
    float amp = frameAmplitude(ibuffer);
    
    if(amp < threshold){
        ++i_low_amp;
    }else{
        i_low_amp = 0;
    }
    
    if(i_frame_session >= session_duration_f && i_low_amp >= minlowampduration_f){
        if(i_session < num_sessions) {
            shifter.setPitchScale(pitch_factors[i_session]);
            i_frame_session=0;
            session_starts[i_session] = i_frame;
            i_session++;
        } else {
            is_finished = i_frame;
        }
    }
    
    memcpy(&input_signal[i_frame*frameSize], ibuffer,sizeof(input_signal[0])*nBufferFrames);
     
    shifter.process(&ibuffer, nBufferFrames,0);
    
    int r = shifter.retrieve(&output_p, nBufferFrames);
    
    memcpy(&output_signal[i_frame*frameSize],output,r*sizeof(float));
    
//     if(add_pink_noise)
//         for (int i=0; i<frameSize; i++ ){  
//             float out = output_signal[i_frame*frameSize+i]+pink_noise[i_noise_frame*frameSize+i];
//             if(out>1.0f) out = 1.0f;
//             if(out<-1.0f) out = -1.0f;
//             *obuffer++ = out;
//             *obuffer++ = out;
//         }
//         ++i_noise_frame;
//         if(i_noise_frame>=num_noise_frames) i_noise_frame=0;
//     else{
//         for (int i=0; i<frameSize; i++ ){       
//             *obuffer++ = output_signal[i_frame*frameSize+i];
//             *obuffer++ = output_signal[i_frame*frameSize+i];
//         }
//     }
        
    for (int i=0; i<frameSize; i++ ){  
            float out = output_signal[i_frame*frameSize+i];
            if(add_pink_noise) out += pink_noise[(i_noise_frame%num_noise_frames)*frameSize+i];
            if(mark_session_starts && i_frame_session < 50) out += (float) 0.05*sin((i_frame*frameSize+i)*(2*3.14159265*1000/sampleRate));
            if(out>1.0f) out = 1.0f;
            if(out<-1.0f) out = -1.0f;
            *obuffer++ = out;
            *obuffer++ = out;
    }
    ++i_noise_frame;
	return 0;
}


void rt_pitch_shifter(void)
{    
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

            mxArray *session_starts_m = mxCreateDoubleMatrix(i_session,1,mxREAL);
            double *session_starts_mp = mxGetPr(session_starts_m);

            for(int i = 0; i<i_session; ++i){
                session_starts_mp[i] = (double) session_starts[i]*frameSize+1;
            }

            if(nlhs >= 1) plhs[0] = input_signal_m;
            if(nlhs >= 2) plhs[1] = output_signal_m;
            if(nlhs >= 3) plhs[2] = session_starts_m;
            
            mexPrintf("Number of processed frames: %i\n",i_frame+1);

            return;
        } else if(v1==0){
            int v2 = is_finished;
            if(v2==0){
                plhs[0] = mxCreateDoubleScalar(0.0);
            } else {
                plhs[0] = mxCreateDoubleScalar((double) (v2*frameSize));
            }
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
            
            if( nrhs < 2 || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetM(prhs[0])>1 && mxGetN(prhs[0])>1)) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notVector", "Second input must be a vector.");
                return;
            }
            
            session_duration_f = v1;
            
            num_sessions = mxGetNumberOfElements(prhs[1]);
            
            double *prhs1 = mxGetPr(prhs[1]);
            
            for(int i = 0; i < num_sessions; ++i){
                pitch_factors[i] = prhs1[i];
            }
            
            if(nrhs>=3){
                if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2])!=1) {
                    mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Third input must be a scalar.");
                    return;
                }
                minlowampduration_f = (int) floor(mxGetScalar(prhs[2]));
            }
            
            if(nrhs>=4){
                if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetNumberOfElements(prhs[3])!=1) {
                    mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Fourth input must be a scalar.");
                    return;
                }
                threshold = (float) mxGetScalar(prhs[3]);
            }
            
            if(nrhs>=5){
                if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxGetNumberOfElements(prhs[4])!=1) {
                    mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Fith input must be a scalar.");
                    return;
                }
                add_pink_noise = (bool) mxGetScalar(prhs[4]);
            }else{
                add_pink_noise = false;
            }
            
            if(nrhs>=6){
                if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || mxGetNumberOfElements(prhs[5])!=1) {
                    mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Fith input must be a scalar.");
                    return;
                }
                noise_gain = (float) mxGetScalar(prhs[5]);
            }
            
            //get noise
            if(add_pink_noise){
                mxArray *num_samples_m = mxCreateDoubleScalar((double) num_noise_frames*frameSize);
                mxArray *pink_noise_m;
                if(mexCallMATLAB(1, &pink_noise_m,1,&num_samples_m,"get_pink_noise")){
                    mexErrMsgIdAndTxt("rt_pitch_shifter:matlabcall", "Couldn't call get_pink_noise function.");
                    return;
                }
                
                if(mxGetNumberOfElements(pink_noise_m)!= num_noise_frames*frameSize){
                    mexErrMsgIdAndTxt("rt_pitch_shifter:matlabcall", "Didn't get enough noise samples.");
                    return;
                }

                //float *pink_noise_mp = (float*) mxGetData(pink_noise_m);
                double *pink_noise_mp = mxGetPr(pink_noise_m);
                for(int i=0; i < num_noise_frames*frameSize; ++i){
                    pink_noise[i] = (float)(noise_gain * pink_noise_mp[i]);
                    //if(pink_noise[i] > 1.0f) pink_noise[i] = 1.0f;
                    //if(pink_noise[i] < -1.0f) pink_noise[i] = -1.0f;
                }
                mxDestroyArray(pink_noise_m);
                mxDestroyArray(num_samples_m);
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
        sum += fabs(samples[i]);
    }
    return sum/frameSize;
        
}