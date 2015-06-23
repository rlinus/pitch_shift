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

double frameAmplitude(double *samples);
double mov_avg(double e);
void gen_beep_sound(double f, double amp);
void gen_piano_sound(double f, double amp);
void gen_drum_sound(double f, double amp);
double draw_var(bool init = false);
void lent_process(double *input, double *output);


const int frameSize = 64;
const int sampleRate = 44100;

//BeeThree instrument2 = BeeThree();

RtAudio dac;

//PitShift lentshifter = PitShift();

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

_dywapitchtracker pitchtracker;

double output[frameSize];
double *output_p = output;
int i_frame;

const size_t data_array_length = 60*44100;
double input_signal[data_array_length];
double output_signal[data_array_length];

double estimated_pitch[data_array_length/frameSize];

double static_factor_sqs[data_array_length/frameSize];
double var_factor_sqs[data_array_length/frameSize];
double control_factor_sqs[data_array_length/frameSize];

const int mov_avg_width = 10;

const double threshold = 0.03;
int voice_onset_f;
int shift_onset_f;
int shift_duration_f;
double pitch_factor;

double control_factor;
double var_factor;
double static_factor;


bool shift_after_voice_onset = false;

int voc_duration_f = 1380;

bool play_ref_sound = true;
bool ref_sound_always_on = false;
bool mark_session_starts = true;
bool mark_session_ends = true;
int start_marker_duration_f = 0.5*sampleRate/(double)frameSize;
int start_marker_onset_f = 30;
int end_marker_duration_f = 1*sampleRate/(double)frameSize;
int ref_sound_duration_f = 0.8*sampleRate/(double)frameSize;

double ref_sound_freq = 125.0;

bool add_pink_noise = false;
double noise_gain = 0.002f;
const int num_noise_frames = 5000;
double pink_noise[num_noise_frames*frameSize];
int i_noise_frame = 0;

const int beep_sound_length = 5*sampleRate;
double beep_sound[beep_sound_length];

const int drum_sound_length = 1*sampleRate;
double drum_sound[drum_sound_length];

const int piano_sound_length = 5*sampleRate;
double piano_sound[piano_sound_length];

volatile int is_finished = 0;

bool do_var = true;
bool do_var_whole_session = false;
double quant_size = 10;
int T_var_f = 8;
double std_dev = 200.0;
double fc = 0.005; //normalized: 1 corresponds to samplingfreq (must be smaller than 0.5)
default_random_engine generator;
normal_distribution<double> distribution_pitch_var(0.0,std_dev);

bool do_control = false;
double kp = 0;
double ki = 0.5;
int control_delay_f = 0.1*sampleRate/(double)frameSize;

double control_error = 0;
double control_error_sum = 0;

void init(void){
    Stk::setSampleRate(sampleRate);
    Stk::setRawwavePath(RAWWAVEPATH);
    
    i_frame = -1;
    
    is_finished = 0;
    voice_onset_f = -1;
    
    control_error_sum = 0;
    //i_noise_frame = 0;
    
    if(shift_after_voice_onset){
        static_factor = 1.0;
    } else {
        static_factor = pitch_factor;
    }
    var_factor = 1;
    control_factor =1;
    
    memset(input_signal,0,sizeof(input_signal));
    memset(output_signal,0,sizeof(output_signal));
    
    draw_var(true);
          
    gen_beep_sound(600,0.1);
    gen_piano_sound(ref_sound_freq, 0.1);
    gen_drum_sound(92.5, 0.1); //92.5: Closed HiHat; 65.4: Base Drum 1

    //init moving average
    for(int i = 0; i<mov_avg_width; ++i){
        mov_avg(0.0f);
    }
    
    dywapitch_inittracking(&pitchtracker);
    
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

	double *ibuffer = (double *) inputBuffer;
	double *obuffer = (double *) outputBuffer;
    
    
    double amp = mov_avg(frameAmplitude(ibuffer));
    
    //ibuffer = &beep_sound[i_frame*frameSize];
    
    memcpy(&input_signal[i_frame*frameSize], ibuffer,sizeof(input_signal[0])*nBufferFrames);
    
    int window_length_factor=16;
    if(i_frame>=window_length_factor){
        estimated_pitch[i_frame] = dywapitch_computepitch(&pitchtracker, &input_signal[(i_frame+1-window_length_factor)*frameSize], 0, window_length_factor*frameSize);
        if (2*estimated_pitch[i_frame]>0.75*ref_sound_freq && 2*estimated_pitch[i_frame]<1.25*ref_sound_freq){
            estimated_pitch[i_frame] = 2 * estimated_pitch[i_frame];
        }else if(0.5*estimated_pitch[i_frame]>0.75*ref_sound_freq && 0.5*estimated_pitch[i_frame]<1.25*ref_sound_freq){
            estimated_pitch[i_frame] = 0.5 * estimated_pitch[i_frame];
        }
    }else{
        estimated_pitch[i_frame] = 0;
    }
    
    if(shift_after_voice_onset){
 
        if(!do_var){
            if(voice_onset_f > -1 && i_frame - voice_onset_f == shift_onset_f){
                static_factor = pitch_factor;
            }
            if(voice_onset_f > -1 && i_frame - voice_onset_f == shift_onset_f + shift_duration_f){
                static_factor = 1.0;
            }
        }else if(!do_var_whole_session){
            static int i_frame_start;
            if(voice_onset_f > -1 && i_frame - voice_onset_f == shift_onset_f){
                i_frame_start = i_frame;
                static_factor = pitch_factor;
            }
            if(voice_onset_f > -1 && i_frame - voice_onset_f >= shift_onset_f && i_frame - voice_onset_f < shift_onset_f + shift_duration_f){
                if((i_frame-i_frame_start)%T_var_f == 0){
                    var_factor = draw_var();
                }
            }
            if(voice_onset_f > -1 && i_frame - voice_onset_f == shift_onset_f + shift_duration_f){
                static_factor = 1.0;
                var_factor = 1.0;
            }
            
        } else{
            if(i_frame%T_var_f == 0){
                var_factor = draw_var();
            }
            if(voice_onset_f > -1 && i_frame - voice_onset_f == shift_onset_f){
                static_factor = pitch_factor;
            }
            if(voice_onset_f > -1 && i_frame - voice_onset_f == shift_onset_f + shift_duration_f){
                static_factor = 1.0;
            }
        }
        
        
    }else if(do_var){
        if(i_frame%T_var_f == 0){
            static_factor = pitch_factor;
            var_factor = draw_var();
        }
    }
    
    if(do_control && voice_onset_f!=-1 && i_frame > voice_onset_f+control_delay_f && (!shift_after_voice_onset || i_frame - voice_onset_f < shift_onset_f + shift_duration_f)){
        if(estimated_pitch[i_frame]>0){
            control_error = 1200*log2(estimated_pitch[i_frame]/ref_sound_freq*static_factor);
        } else{
            control_error = 0;
        }
        control_error_sum += control_error*frameSize/(double)sampleRate;
        double control_factor_cents = kp*control_error+ki*control_error_sum;
        
        const double max_control_error = 200;
        if(control_factor_cents > max_control_error) control_factor_cents = max_control_error;
        if(control_factor_cents < -1*max_control_error) control_factor_cents = -1*max_control_error;
        control_factor = pow(2,control_factor_cents/1200);
    }else{
        control_factor = 1;
    }
    
    if (voice_onset_f == -1 && amp > threshold && i_frame > ref_sound_duration_f){
            voice_onset_f = i_frame;
    }
    
    if(voice_onset_f > -1 && i_frame == voice_onset_f+voc_duration_f+end_marker_duration_f*mark_session_ends){
        is_finished = 1;
    }
    
    static_factor_sqs[i_frame] = static_factor;
    var_factor_sqs[i_frame] = var_factor;
    control_factor_sqs[i_frame] = control_factor;
    
    #if defined CPVSHIFTER
    	cpvPitchShift(static_factor*var_factor*control_factor, ibuffer, output);
    #elif defined SMBSHIFTER
        smbPitchShift(static_factor*var_factor*control_factor, ibuffer, output);
    #elif defined DIRACSHIFTER
        float in[frameSize];
        float out[frameSize];
        for(int i=0; i < frameSize; ++i){
            in[i] = (float) ibuffer[i];
        }
        
        DiracFxProcessFloatInterleaved(1.0, static_factor*var_factor*control_factor, in, out, frameSize, diracFx);
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
        
        shifter.setPitchScale(static_factor*var_factor*control_factor);
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
    
    //if(i_frame>=window_length_factor*4){
    memcpy(&output_signal[i_frame*frameSize],output,frameSize*sizeof(double));
    //}
    
     //memcpy(&output_signal[i_frame*frameSize],ibuffer,frameSize*sizeof(double));
    
        
    for (int i=0; i<frameSize; i++ ){  
            double out = output_signal[i_frame*frameSize+i];
            if(add_pink_noise) out += pink_noise[(i_noise_frame%num_noise_frames)*frameSize+i];
            if(play_ref_sound && (i_frame < ref_sound_duration_f || ref_sound_always_on)){
                out += (double) piano_sound[((i_frame)*frameSize+i)%piano_sound_length];
            }
            if(mark_session_starts && i_frame >= start_marker_onset_f && i_frame < start_marker_onset_f+start_marker_duration_f){
                out += (double) drum_sound[(i_frame-start_marker_onset_f)*frameSize+i];
            }
            if(mark_session_ends && voice_onset_f > -1 && i_frame >= voice_onset_f+voc_duration_f && i_frame < voice_onset_f+voc_duration_f+end_marker_duration_f){
                out += (double) drum_sound[(i_frame-(voice_onset_f+voc_duration_f))*frameSize+i];
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
    oparameters.deviceId = DEVICEID; //dac.getDefaultOutputDevice();
	oparameters.nChannels = 2;

	RtAudio::StreamParameters iparameters;
	iparameters.deviceId = DEVICEID; //dac.getDefaultInputDevice();
	iparameters.nChannels = 1;
    
    RtAudio::StreamOptions options;
    options.flags = RTAUDIO_HOG_DEVICE;// | RTAUDIO_SCHEDULE_REALTIME | RTAUDIO_MINIMIZE_LATENCY;
    //options.priority = 10000000;
    
    unsigned int bufferFrames = frameSize;
    
    
	try {
		dac.openStream( &oparameters, &iparameters, RTAUDIO_FLOAT64, (unsigned int)Stk::sampleRate(), &bufferFrames, &tick, NULL, &options);
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

            mxArray *input_signal_m = mxCreateDoubleMatrix((i_frame+1)*frameSize,1,mxREAL);
            double *input_signal_mp = mxGetPr(input_signal_m);
            mxArray *output_signal_m = mxCreateDoubleMatrix((i_frame+1)*frameSize,1,mxREAL);
            double *output_signal_mp = mxGetPr(output_signal_m);

            for(int i = 0; i<(i_frame+1)*frameSize; ++i){
                input_signal_mp[i] = (double) input_signal[i];
                output_signal_mp[i] = (double) output_signal[i];
            }

            mxArray *voice_onset_m = mxCreateDoubleScalar(voice_onset_f);
            
            mxArray *static_factor_sqs_m = mxCreateDoubleMatrix(i_frame+1,1,mxREAL);
            double *static_factor_sqs_mp = mxGetPr(static_factor_sqs_m);
            mxArray *var_factor_sqs_m = mxCreateDoubleMatrix(i_frame+1,1,mxREAL);
            double *var_factor_sqs_mp = mxGetPr(var_factor_sqs_m);
            mxArray *control_factor_sqs_m = mxCreateDoubleMatrix(i_frame+1,1,mxREAL);
            double *control_factor_sqs_mp = mxGetPr(control_factor_sqs_m);
            mxArray *estimated_pitch_m = mxCreateDoubleMatrix(i_frame+1,1,mxREAL);
            double *estimated_pitch_mp = mxGetPr(estimated_pitch_m);
            
            for(int i = 0; i<i_frame+1; ++i){
                static_factor_sqs_mp[i] = (double) static_factor_sqs[i];
                control_factor_sqs_mp[i] = (double) control_factor_sqs[i];
                var_factor_sqs_mp[i] = (double) var_factor_sqs[i];
                estimated_pitch_mp[i] = (double) estimated_pitch[i];
            }

            if(nlhs >= 1) plhs[0] = input_signal_m;
            if(nlhs >= 2) plhs[1] = output_signal_m;
            if(nlhs >= 3) plhs[2] = voice_onset_m;
            if(nlhs >= 4) plhs[3] = static_factor_sqs_m;
            if(nlhs >= 5) plhs[4] = var_factor_sqs_m;
            if(nlhs >= 6) plhs[5] = control_factor_sqs_m;
            if(nlhs >= 7) plhs[6] = estimated_pitch_m;
            
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
            
            if( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) || mxGetNumberOfElements(prhs[6])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Seventh input must be a scalar.");
                return;
            }
            do_var = (int) floor(mxGetScalar(prhs[6]));
            
            if( !mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) || mxGetNumberOfElements(prhs[7])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Eighth input must be a scalar.");
                return;
            }
            std_dev = (double) mxGetScalar(prhs[7]);
            
            if( !mxIsDouble(prhs[8]) || mxIsComplex(prhs[8]) || mxGetNumberOfElements(prhs[8])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Nineth input must be a scalar.");
                return;
            }
            fc = (double) mxGetScalar(prhs[8]);
            
            if( !mxIsDouble(prhs[9]) || mxIsComplex(prhs[9]) || mxGetNumberOfElements(prhs[9])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "Tenth input must be a scalar.");
                return;
            }
            do_control = (int) floor(mxGetScalar(prhs[9]));
            
            if( !mxIsDouble(prhs[10]) || mxIsComplex(prhs[10]) || mxGetNumberOfElements(prhs[10])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "11th input must be a scalar.");
                return;
            }
            kp = (double) mxGetScalar(prhs[10]);
            
            if( !mxIsDouble(prhs[11]) || mxIsComplex(prhs[11]) || mxGetNumberOfElements(prhs[11])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "12th input must be a scalar.");
                return;
            }
            ki = (double) mxGetScalar(prhs[11]);
            
            if( !mxIsDouble(prhs[12]) || mxIsComplex(prhs[12]) || mxGetNumberOfElements(prhs[12])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "13th input must be a scalar.");
                return;
            }
            do_var_whole_session = (int) mxGetScalar(prhs[12]);
            
            if( !mxIsDouble(prhs[13]) || mxIsComplex(prhs[13]) || mxGetNumberOfElements(prhs[13])!=1) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:notScalar", "14th input must be a scalar.");
                return;
            }
            noise_gain = (double) mxGetScalar(prhs[13]);
            if(noise_gain > 0){
                add_pink_noise = true;
            }else{
                add_pink_noise = false;
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

                //double *pink_noise_mp = (double*) mxGetData(pink_noise_m);
                double *pink_noise_mp = mxGetPr(pink_noise_m);
                for(int i=0; i < num_noise_frames*frameSize; ++i){
                    pink_noise[i] = (double)(noise_gain * pink_noise_mp[i]);
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

double frameAmplitude(double *samples) {
    double sum = 0;
    for(int i = 0; i < frameSize; ++i) {
        sum += fabs(samples[i]);
    }
    return sum/frameSize;
        
}


double mov_avg(double e){
    static double b[mov_avg_width];
    static int p = 0;
    
    b[p] = e;
    p++;
    if(p >= mov_avg_width) p = 0;
    
    double sum = 0;
    for(int i = 0; i < mov_avg_width; ++i){
        sum += b[i];
    }
    return sum/mov_avg_width;
}

void gen_beep_sound(double f, double amp){
    const double pi = 3.14159265358979323846;
    for(int i = 0; i < beep_sound_length; ++i){
        beep_sound[i] = amp*sin(i *f*2*pi/sampleRate);
    }
    return;
}

void gen_piano_sound(double f, double amp){
    try{
        Instrmnt *instrument = new Rhodey();
        //Instrmnt *instrument = new BeeThree();  //causes unknown exception sometimes (but never after compiling
        //Instrmnt *instrument = new Wurley();
        instrument->noteOn(f, amp);

        for(int i = 0; i < piano_sound_length; ++i){
            piano_sound[i] = (double) instrument->tick();
        }

        delete instrument;
    }catch(...){
        mexPrintf("Problem generating piano sound. Recompile!\n");
    }
}

void gen_drum_sound(double f, double amp){
    try{
        Instrmnt *instrument = new Drummer();
        instrument->noteOn(f, amp);

        for(int i = 0; i < drum_sound_length; ++i){
            drum_sound[i] = (double) instrument->tick();
        }

        delete instrument;
    }catch(...){
        mexPrintf("Problem generating drum sound. Recompile!\n");
    }
}

// double draw_var(bool init =false){
//     static double v = 0;
//     v = (1-delta) * v + delta * distribution_pitch_var(generator);
//     //double v = distribution_pitch_var(generator);
//     //if(v > 2*std_var) v=2*std_var;
//     //if(v < -2*std_var) v=-2*std_var;
//     return pow(2.0,v/1200.0);
// }

double draw_var(bool init){
    //second order butterworth filter (see: http://www.kwon3d.com/theory/filtering/fil.html)
    
    static double ym1 = 0;
    static double ym2 = 0;
    static double xm1 = 0;
    static double xm2 = 0;
    
    const double pi = 3.14159265358979323846;
    
    static double omega;
    static double c;
    
    static double a0;
    static double a1;
    static double a2;
    static double b1;
    static double b2;
    
    
    if(init){
        ym1 = 0;
        ym2 = 0;
        xm1 = 0;
        xm2 = 0;
        
        omega = tan(pi*fc);
        c = 1+sqrt(2)*omega+omega*omega;

        a0 = omega*omega/c;
        a1 = 2*a0;
        a2 = a0;
        b1 = 2*(omega*omega-1)/c;
        b2 = (1-sqrt(2)*omega+omega*omega)/c;
        
        return 0;
    }
    
    double x = distribution_pitch_var(generator);
    
    double y;
    
    if(fc>=0.5){
        y=x;
    }else{
        y = a0*x+a1*xm1+a2*xm2-b1*ym1-b2*ym2;

        ym2 = ym1;
        ym1 = y;
        xm2 = xm1;
        xm1 = x;
    }
    
    if(quant_size > 0){
        y = quant_size*round(y/quant_size);
    }
    
    return pow(2.0,y/1200.0);
    
}

// void lent_process(double *input, double *output){
//     for(int i = 0; i<frameSize; ++i){
//         output[i] = (double) lentshifter.tick((StkFloat) input[i]);
//     }
//     return;
// }
