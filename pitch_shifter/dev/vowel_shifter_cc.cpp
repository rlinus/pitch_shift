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

double frameAmplitude(float *samples);
double mov_avg(double e);
void gen_sin(float *p,double f, double amp, int n);

const int frameSize = 64;
const int sampleRate = 44100; 

RtAudio dac;

RubberBandStretcher::Options options =  RubberBandStretcher::OptionProcessRealTime |
                                        RubberBandStretcher::OptionFormantPreserved |
                                        RubberBandStretcher::OptionWindowShort |
                                        RubberBandStretcher::OptionTransientsSmooth |
                                        RubberBandStretcher::OptionPitchHighConsistency;
RubberBandStretcher shifter(sampleRate,1,options,1.0,1.0);

float output[frameSize];
float *output_p = output;
int i_frame;

const size_t data_array_length = 60*44100;
float input_signal[data_array_length];
float output_signal[data_array_length];

const int mov_avg_width = 5;

const double threshold = 0.005f;
int voice_onset_f;
int shift_onset_f;
int shift_duration_f;
double pitch_factor;


int voc_duration_f = 1380;
int marker_duration_f = 100;
bool mark_session_starts = true;
bool mark_session_ends = true;

bool add_pink_noise = false;
float noise_gain = 0.002f;
const int num_noise_frames = 5000;
float pink_noise[num_noise_frames*frameSize];
int i_noise_frame = 0;

float start_marker[sampleRate];

volatile int is_finished = 0;

void init_global_vars(void){
    i_frame = -1;
    
    is_finished = 0;
    voice_onset_f = -1;
    
    //i_noise_frame = 0;
    
    shifter.reset();
    shifter.setPitchScale(pitch_factor);
    
    memset(input_signal,0,sizeof(input_signal));
    memset(output_signal,0,sizeof(output_signal));
    
    gen_sin(start_marker,500,0.03,sampleRate);
    
    for(int i = 0; i<mov_avg_width; ++i){
        mov_avg(0.0f);
    }
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

	float *ibuffer = (float *) inputBuffer;
	float *obuffer = (float *) outputBuffer;
    
    
    double amp = mov_avg(frameAmplitude(ibuffer));
    
//     if (voice_onset_f == -1 && i_frame == 80){
//         shifter.setPitchScale(1.0);
//     }
//     
    if (voice_onset_f == -1 && amp > threshold){
        voice_onset_f = i_frame;
        //shifter.setPitchScale(1.0);
    }
    
//     if(voice_onset_f > -1 && i_frame - voice_onset_f == shift_onset_f){
//         shifter.setPitchScale(pitch_factor);
//         //mexPrintf("pitch: %f,frame: %i\n",amp, i_frame);
//     }
   
//     if(voice_onset_f > -1 && i_frame - voice_onset_f == shift_onset_f + shift_duration_f){
//         shifter.setPitchScale(1.0);
//         //mexPrintf("pitch: %f,frame: %i\n",amp, i_frame);
//     }
    
    if(voice_onset_f > -1 && i_frame == voice_onset_f+voc_duration_f+marker_duration_f){
        is_finished = 1;
    }
    
    memcpy(&input_signal[i_frame*frameSize], ibuffer,sizeof(input_signal[0])*nBufferFrames);
     
    shifter.process(&ibuffer, nBufferFrames,0);
    
    int r = shifter.retrieve(&output_p, nBufferFrames);
    
    //memcpy(&output_signal[i_frame*frameSize],output,r*sizeof(float));
    
    if(voice_onset_f > -1 && i_frame - voice_onset_f >= shift_onset_f && i_frame - voice_onset_f <= shift_onset_f + shift_duration_f){
        memcpy(&output_signal[i_frame*frameSize],output,r*sizeof(float));
    }else{
        memcpy(&output_signal[i_frame*frameSize],ibuffer,sizeof(input_signal[0])*nBufferFrames);
    }
        
    for (int i=0; i<frameSize; i++ ){  
            float out = output_signal[i_frame*frameSize+i];
            if(add_pink_noise) out += pink_noise[(i_noise_frame%num_noise_frames)*frameSize+i];
            if(mark_session_starts && i_frame >= 100 && i_frame < 100+marker_duration_f){
                out += (float) start_marker[(i_frame-100)*frameSize+i];
            }
            else if(mark_session_ends && voice_onset_f > -1 && i_frame >= voice_onset_f+voc_duration_f && i_frame < voice_onset_f+voc_duration_f+marker_duration_f){
                out += (float) start_marker[(i_frame-(voice_onset_f+voc_duration_f))*frameSize+i];
            }
            if(out>1.0f) out = 1.0f;
            if(out<-1.0f) out = -1.0f;
            *obuffer++ = out;
            *obuffer++ = out;
    }
    ++i_noise_frame;
	return 0;
}

void start_stream(void)
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

            mxArray *voice_onset_m = mxCreateDoubleScalar(voice_onset_f);

            if(nlhs >= 1) plhs[0] = input_signal_m;
            if(nlhs >= 2) plhs[1] = output_signal_m;
            if(nlhs >= 3) plhs[2] = voice_onset_m;
            
            mexPrintf("Number of processed frames: %i\n",i_frame+1);

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
            
            if( nrhs < 4) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:missinginput", "Need at least two input arguments.");
                return;
            }

            if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfElements(prhs[0])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "First input must be a scalar.");
                return;
            }
            shift_onset_f = (int) floor(mxGetScalar(prhs[0]));

            if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetNumberOfElements(prhs[1])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Second input must be a scalar.");
                return;
            }
            pitch_factor = mxGetScalar(prhs[1]);
            
            if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Third input must be a scalar.");
                return;
            }
            shift_duration_f = (int) floor(mxGetScalar(prhs[2]));
            
            if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetNumberOfElements(prhs[3])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Third input must be a scalar.");
                return;
            }
            voc_duration_f = (int) floor(mxGetScalar(prhs[3]));
            
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
            start_stream();
        }
    }
    
    return;
}

double frameAmplitude(float *samples) {
    double sum = 0;
    for(int i = 0; i < frameSize; ++i) {
        sum += fabs(samples[i]);
    }
    return sum/frameSize;
        
}


double mov_avg(double e){
    static float b[mov_avg_width];
    static int p = 0;
    
    b[p] = e;
    p++;
    if(p >= mov_avg_width) p = 0;
    
    float sum = 0;
    for(int i = 0; i < mov_avg_width; ++i){
        sum += b[i];
    }
    return sum/mov_avg_width;
}

void gen_sin(float *p,double f, double amp, int n){
    const double pi = 3.14159265358979323846;
    for(int i = 1; i < n; ++i){
        p[i] =(float) amp*sin(i *f*2*pi/sampleRate);
    }
    return;
}