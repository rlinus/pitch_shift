#include <iostream>
#include <stdlib.h>
#include <random>
#include <string>
#include "mex.h"
#include "Stk.h"
#include "RtAudio.h"
#include "RubberBandStretcher.h"
#include "Wurley.h"
#include "Rhodey.h"
#include "BeeThree.h"
#include "Drummer.h"
//#include "PitShift.h"
//#include "LentPitShift.h"

using namespace std;
using namespace stk;
using namespace RubberBand;

double frameAmplitude(float *samples);
double mov_avg(double e);
void gen_beep_sound(double f, double amp);
void gen_piano_sound(double f, double amp);
void gen_drum_sound(double f, double amp);
double draw_var(void);
void lent_process(float *input, float *output);

const int frameSize = 64;
const int sampleRate = 44100;

//BeeThree instrument2 = BeeThree();

RtAudio dac;

RubberBandStretcher::Options options =  RubberBandStretcher::OptionProcessRealTime |
                                        RubberBandStretcher::OptionFormantPreserved |
                                        RubberBandStretcher::OptionWindowShort |
                                        //RubberBandStretcher::OptionSmoothingOn |
                                        //RubberBandStretcher::OptionPhaseIndependent |
                                        RubberBandStretcher::OptionTransientsSmooth |
                                        RubberBandStretcher::OptionPitchHighQuality;//Consistency;
RubberBandStretcher shifter(sampleRate,1,options,1.0,1.0);

//PitShift lentshifter = PitShift();

float output[frameSize];
float *output_p = output;
int i_frame;

const size_t data_array_length = 60*44100;
float input_signal[data_array_length];
float output_signal[data_array_length];

int num_sample_available[data_array_length/frameSize];
int sw_latency[data_array_length/frameSize];

const int mov_avg_width = 5;

const double threshold = 0.02f;
int voice_onset_f;
int shift_onset_f;
int shift_duration_f;
double pitch_factor;

bool shift_after_voice_onset = false;

int voc_duration_f = 1380;

bool play_ref_sound = true;
bool ref_sound_always_on = false;
bool mark_session_starts = true;
bool mark_session_ends = true;
int start_marker_duration_f = 350;
int start_marker_onset_f = 30;
int end_marker_duration_f = 350;
int ref_sound_duration_f = 350;

double ref_sound_freq = 125.0;

bool add_pink_noise = false;
float noise_gain = 0.002f;
const int num_noise_frames = 5000;
float pink_noise[num_noise_frames*frameSize];
int i_noise_frame = 0;

const int beep_sound_length = 1*sampleRate;
float beep_sound[beep_sound_length];

const int drum_sound_length = 1*sampleRate;
float drum_sound[drum_sound_length];

const int piano_sound_length = 5*sampleRate;
float piano_sound[piano_sound_length];

volatile int is_finished = 0;

bool do_var = false;
double std_dev = 50.0;
double delta = 1;
int T_var_f = 64;
default_random_engine generator;
normal_distribution<double> distribution_pitch_var(0.0,std_dev);

void init(void){
    Stk::setSampleRate(sampleRate);
    Stk::setRawwavePath( "../stk-4.5.0/rawwaves/" );
    
    i_frame = -1;
    
    is_finished = 0;
    voice_onset_f = -1;
    
    //i_noise_frame = 0;
    
    shifter.reset();
    if(shift_after_voice_onset){
        shifter.setPitchScale(1.0);
        //lentshifter.setShift(1.0);
    } else {
        shifter.setPitchScale(pitch_factor);
        //lentshifter.setShift(pitch_factor);
    }
    
    memset(input_signal,0,sizeof(input_signal));
    memset(output_signal,0,sizeof(output_signal));
          
    
    gen_beep_sound(600,0.1);
    gen_piano_sound(ref_sound_freq, 0.1);
    gen_drum_sound(92.5, 0.25); //92.5: Closed HiHat; 65.4: Base Drum 1

    //init moving average
    for(int i = 0; i<mov_avg_width; ++i){
        mov_avg(0.0f);
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

	float *ibuffer = (float *) inputBuffer;
	float *obuffer = (float *) outputBuffer;
    
    
    double amp = mov_avg(frameAmplitude(ibuffer));
    
    if(shift_after_voice_onset){
 
        if(!do_var){
            if(voice_onset_f > -1 && i_frame - voice_onset_f == shift_onset_f){
                shifter.setPitchScale(pitch_factor);
                //lentshifter.setShift(pitch_factor);
                //mexPrintf("pitch: %f,frame: %i\n",amp, i_frame);
            } 
        }else{
            static int i_frame_start;
            if(voice_onset_f > -1 && i_frame - voice_onset_f == shift_onset_f){
                i_frame_start = i_frame;
            }
            if(voice_onset_f > -1 && i_frame - voice_onset_f >= shift_onset_f && i_frame - voice_onset_f < shift_onset_f + shift_duration_f){
                if((i_frame-i_frame_start)%T_var_f == 0){
                    shifter.setPitchScale(pitch_factor*draw_var());
                    //lentshifter.setShift(pitch_factor*draw_var());
                }
            }
            
        }
        
        if(voice_onset_f > -1 && i_frame - voice_onset_f == shift_onset_f + shift_duration_f){
                shifter.setPitchScale(1.0);
                //lentshifter.setShift(1.0);
                //mexPrintf("pitch: %f,frame: %i\n",amp, i_frame);
        }
    }else if(do_var){
        if(i_frame%T_var_f == 0){
            shifter.setPitchScale(pitch_factor*draw_var());
            //lentshifter.setShift(pitch_factor*draw_var());
        }
    }
    
    if (voice_onset_f == -1 && amp > threshold){
            voice_onset_f = i_frame;
    }
    
    if(voice_onset_f > -1 && i_frame == voice_onset_f+voc_duration_f+end_marker_duration_f){
        is_finished = 1;
    }
    
    memcpy(&input_signal[i_frame*frameSize], ibuffer,sizeof(input_signal[0])*nBufferFrames);
     
    shifter.process(&ibuffer, nBufferFrames,0);
    
    int r = shifter.available();
    num_sample_available[i_frame] = r;
    sw_latency[i_frame] = shifter.getLatency();
    
    //the latency of the shifter is 993 samples for pitch_factor = 1.0 (measured). I add two frames latency, to get consitency after shifts, therefore the software latency is 993+2*64=1121 samples
    if(r>=frameSize && i_frame+1>=18){
        shifter.retrieve(&output_p, frameSize);
        memcpy(&output_signal[i_frame*frameSize],output,frameSize*sizeof(float));
    }
    
//     lent_process(ibuffer,output_p);
//     memcpy(&output_signal[i_frame*frameSize],output,frameSize*sizeof(float));
    
        
    for (int i=0; i<frameSize; i++ ){  
            float out = output_signal[i_frame*frameSize+i];
            if(add_pink_noise) out += pink_noise[(i_noise_frame%num_noise_frames)*frameSize+i];
            if(play_ref_sound && (i_frame < ref_sound_duration_f || ref_sound_always_on)){
                out += (float) piano_sound[((i_frame)*frameSize+i)%piano_sound_length];
            }
            if(mark_session_starts && i_frame >= start_marker_onset_f && i_frame < start_marker_onset_f+start_marker_duration_f){
                out += (float) drum_sound[(i_frame-start_marker_onset_f)*frameSize+i];
            }
            if(mark_session_ends && voice_onset_f > -1 && i_frame >= voice_onset_f+voc_duration_f && i_frame < voice_onset_f+voc_duration_f+end_marker_duration_f){
                out += (float) drum_sound[(i_frame-(voice_onset_f+voc_duration_f))*frameSize+i];
            }
            //out *= 2;
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
    
    RtAudio::StreamParameters oparameters;
    oparameters.deviceId = 0; //dac.getDefaultOutputDevice();
	oparameters.nChannels = 2;

	RtAudio::StreamParameters iparameters;
	iparameters.deviceId = 0; //dac.getDefaultInputDevice();
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

            mxArray *voice_onset_m = mxCreateDoubleScalar(voice_onset_f);
            
            mxArray *num_sample_available_m = mxCreateDoubleMatrix(i_frame+1,1,mxREAL);
            double *num_sample_available_mp = mxGetPr(num_sample_available_m);
            mxArray *sw_latency_m = mxCreateDoubleMatrix(i_frame+1,1,mxREAL);
            double *sw_latency_mp = mxGetPr(sw_latency_m);
            
            for(int i = 0; i<i_frame+1; ++i){
                num_sample_available_mp[i] = (double) num_sample_available[i];
                sw_latency_mp[i] = (double) sw_latency[i];
            }

            if(nlhs >= 1) plhs[0] = input_signal_m;
            if(nlhs >= 2) plhs[1] = output_signal_m;
            if(nlhs >= 3) plhs[2] = voice_onset_m;
            if(nlhs >= 4) plhs[3] = num_sample_available_m;
            if(nlhs >= 5) plhs[4] = sw_latency_m;
            
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
            
            if(v1==1){
                shift_after_voice_onset = false;
            }else if(v1==2){
                shift_after_voice_onset = true;
            }

            if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetNumberOfElements(prhs[1])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Second input must be a scalar.");
                return;
            }
            pitch_factor = mxGetScalar(prhs[1]);
            
            if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Third input must be a scalar.");
                return;
            }
            voc_duration_f = (int) floor(mxGetScalar(prhs[2]));
            
            if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetNumberOfElements(prhs[3])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Fourth input must be a scalar.");
                return;
            }
            double f = mxGetScalar(prhs[3]);
            if((int)f==0){
                play_ref_sound = false;
                mark_session_starts = true;
                start_marker_onset_f = 30;
            }else{
                play_ref_sound = true;
                if(f<0){
                    ref_sound_always_on = true;
                    mark_session_starts = true;
                    start_marker_onset_f = 670;
                }else{
                    ref_sound_always_on = false;
                    mark_session_starts = false;
                }
            }
            ref_sound_freq = fabs(f);
            
            
            if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxGetNumberOfElements(prhs[4])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Fivth input must be a scalar.");
                return;
            }
            shift_onset_f = (int) floor(mxGetScalar(prhs[4]));
            
            if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || mxGetNumberOfElements(prhs[5])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Sixth input must be a scalar.");
                return;
            }
            shift_duration_f = (int) floor(mxGetScalar(prhs[5]));
            
            
            
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
            
            
            init();
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

void gen_beep_sound(double f, double amp){
    const double pi = 3.14159265358979323846;
    for(int i = 0; i < beep_sound_length; ++i){
        beep_sound[i] =(float) amp*sin(i *f*2*pi/sampleRate);
    }
    return;
}

void gen_piano_sound(double f, double amp){
      
    //Instrmnt *instrument = new Rhodey();
    Instrmnt *instrument = new BeeThree();  //causes unknown exception sometimes (but never after compiling
    //Instrmnt *instrument = new Wurley();
    instrument->noteOn(f, amp);
    
    for(int i = 0; i < piano_sound_length; ++i){
        piano_sound[i] = (float) instrument->tick();
    }
    
    delete instrument;
}

void gen_drum_sound(double f, double amp){
    Instrmnt *instrument = new Drummer();
    instrument->noteOn(f, amp);
    
    for(int i = 0; i < drum_sound_length; ++i){
        drum_sound[i] = (float) instrument->tick();
    }
    
    delete instrument;
}

double draw_var(void){
    static double v = 0;
    v = (1-delta) * v + delta * distribution_pitch_var(generator);
    //double v = distribution_pitch_var(generator);
    //if(v > 2*std_var) v=2*std_var;
    //if(v < -2*std_var) v=-2*std_var;
    return pow(2.0,v/1200.0);
}

// void lent_process(float *input, float *output){
//     for(int i = 0; i<frameSize; ++i){
//         output[i] = (float) lentshifter.tick((StkFloat) input[i]);
//     }
//     return;
// }
