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
#include <minmax.h>



#include "cpvPitchShift.h"
#include "smbPitchShift.h"

#define RUBBERBAND //if RUBBERBAND lib is available

#if defined RUBBERBAND
#include "RubberBandStretcher.h"
using namespace RubberBand;
#endif



using namespace std;
using namespace stk;

int shifterId = 1;
int deviceId = 0;

double frameAmplitude(double *samples);
double mov_avg(double e);
double sig_power(void);
void gen_beep_sound(double f, double amp);
void gen_piano_sound(double f, double amp);
void gen_drum_sound(double f, double amp);
double draw_var(bool init = false);
std::string get_mex_path(void);


const int frameSize = 64;
const int sampleRate = 44100;

//BeeThree instrument2 = BeeThree();

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

_dywapitchtracker pitchtracker;

double output[frameSize];
double *output_p = output;
int i_frame;

const size_t data_array_length = 120*sampleRate;
double input_signal[data_array_length];
double output_signal[data_array_length];

double estimated_pitch[data_array_length/frameSize];

double static_factor_sqs[data_array_length/frameSize];
double var_factor_sqs[data_array_length/frameSize];
double control_factor_sqs[data_array_length/frameSize];

const int mov_avg_width = 100;
const int sig_power_width = 30;


double start_threshold = 0.01;
double stop_threshold = 0.01;

int voice_onset_f;
int shift_onset_f;
int shift_duration_f;
double pitch_factor;

double control_factor;
double var_factor;
double static_factor;


bool shift_after_voice_onset = false;

int voc_duration_f = 1380;
double feedback_gain = 1;

bool play_ref_sound = true;
bool ref_sound_always_on = false;
bool mark_session_starts = false;
bool mark_session_ends = false;
int start_marker_duration_f = 0.5*sampleRate/(double)frameSize;
int start_marker_onset_f = 30;
int end_marker_duration_f = 1*sampleRate/(double)frameSize;
int ref_sound_duration_f = 1*sampleRate/(double)frameSize;

bool custom_ref_sound = false;
double ref_sound_freq = 125.0;
double ref_sound_amplitude = 0.1;;

bool add_pink_noise = false;
double noise_gain = 0.005;
const int num_noise_frames = 5000;
double pink_noise[num_noise_frames*frameSize];

bool adaptive_noise_level = false;
double min_noise_level = 0.005;
double max_noise_level = 0.05;

const int beep_sound_length = 5*sampleRate;
double beep_sound[beep_sound_length];

const int drum_sound_length = 1*sampleRate;
double drum_sound[drum_sound_length];

const int ref_sound_length = 60*sampleRate;
double ref_sound[ref_sound_length];

volatile int is_finished = 0;

bool do_var = true;
bool do_var_whole_session = false;
double var_quant_size = 10;
int T_var_f = 8;
double std_dev = 200.0;
double fc = 0.005; //normalized: 1 corresponds to nyquistfreq (must be smaller than 1)
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
    Stk::setRawwavePath(get_mex_path()+"rawwaves/");
    //mexPrintf("%s\n",(get_mex_path()+"../stk-4.5.0/rawwaves/").c_str());
    
    i_frame = -1;
    
    is_finished = 0;
    voice_onset_f = -1;
    
    control_error_sum = 0;
    
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
          
    if(play_ref_sound && !custom_ref_sound) gen_piano_sound(ref_sound_freq,ref_sound_amplitude);
    //gen_beep_sound(600,0.1);
    //gen_drum_sound(92.5, 0.2); //92.5: Closed HiHat; 65.4: Base Drum 1

    //init moving average
    for(int i = 0; i<mov_avg_width; ++i){
        mov_avg(0.0f);
    }
    
    dywapitch_inittracking(&pitchtracker);
    
    switch(shifterId){
        case 0:
            smbPitchShiftInit(frameSize, 1024/frameSize, sampleRate);
            break;
        case 1:
            cpvPitchShiftInit(frameSize, 1024/frameSize, sampleRate);
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

	double *ibuffer = (double *) inputBuffer;
	double *obuffer = (double *) outputBuffer;
    
    double frame_amp = frameAmplitude(ibuffer);
    double amp = mov_avg(frame_amp);
    
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
    
    if (voice_onset_f == -1 && amp > start_threshold && i_frame > ref_sound_duration_f){
            voice_onset_f = i_frame;
    }
    
    if(voice_onset_f > -1 && i_frame == voice_onset_f+voc_duration_f+end_marker_duration_f*mark_session_ends){
        is_finished = 1;
    }
    if(is_finished == 1 && amp < stop_threshold){
        is_finished = 2;
    }
    
    static_factor_sqs[i_frame] = static_factor;
    var_factor_sqs[i_frame] = var_factor;
    control_factor_sqs[i_frame] = control_factor;

    switch(shifterId){
        case 0:
            smbPitchShift(static_factor*var_factor*control_factor, ibuffer, output);
            break;
        case 1:
            cpvPitchShift(static_factor*var_factor*control_factor, ibuffer, output);
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
            break;
        #endif
    }

    
    //if(i_frame>=window_length_factor*4){
    memcpy(&output_signal[i_frame*frameSize],output,frameSize*sizeof(double));
    //}
    
     //memcpy(&output_signal[i_frame*frameSize],ibuffer,frameSize*sizeof(double));
    
        
    for (int i=0; i<frameSize; i++ ){  
            double out = feedback_gain*output_signal[i_frame*frameSize+i];
            if(add_pink_noise && (!play_ref_sound || i_frame > ref_sound_duration_f)){
                if(adaptive_noise_level){
                    out += pink_noise[(i_frame%num_noise_frames)*frameSize+i]*min(max(sig_power(),min_noise_level),max_noise_level);
                }else{
                    out += pink_noise[(i_frame%num_noise_frames)*frameSize+i];
                }
            }
            if(play_ref_sound && (i_frame < ref_sound_duration_f || ref_sound_always_on)){
                out += (double) ref_sound[((i_frame)*frameSize+i)%ref_sound_length];
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
	return 0;
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
    options.flags = RTAUDIO_HOG_DEVICE | RTAUDIO_MINIMIZE_LATENCY;// | RTAUDIO_SCHEDULE_REALTIME;
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

            mxArray *voice_onset_m = mxCreateDoubleScalar(voice_onset_f*frameSize/double(sampleRate));
            
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
            
            if( nrhs < 2) {
                mexErrMsgIdAndTxt("rt_pitch_shifter:missinginput", "Need at least two input arguments.");
                return;
            }
            
            if(!mxIsStruct(prhs[1])){
                mexErrMsgIdAndTxt("rt_pitch_shifter:wronginput", "Second argument must be a struct.");
                return;
            }
            
            //get input parameters
            mxArray * fieldptr;
            
            fieldptr = mxGetField(prhs[1], 0, "shift_full_trial");
            if(fieldptr){
                shift_after_voice_onset = !(bool)mxGetScalar(fieldptr);
            }else {
                shift_after_voice_onset = 1;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "pitch_factor");
            if(fieldptr){
                pitch_factor = mxGetScalar(fieldptr);
            }else{
                pitch_factor = 1;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "voc_duration");
            if(fieldptr){
                voc_duration_f = round(mxGetScalar(fieldptr)*double(sampleRate)/frameSize);
            }else{
                voc_duration_f = sampleRate/frameSize;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "shift_onset");
            if(fieldptr){
                shift_onset_f = round(mxGetScalar(fieldptr)*double(sampleRate)/frameSize);
            }else{
                shift_onset_f = voc_duration_f/4;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "shift_duration");
            if(fieldptr){
                shift_duration_f = round(mxGetScalar(fieldptr)*double(sampleRate)/frameSize);
            }else{
                shift_duration_f = voc_duration_f/4;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "play_ref_sound");
            if(fieldptr){
                play_ref_sound = (bool) mxGetScalar(fieldptr);
            }else{
                play_ref_sound = 0;
            }
            
            
            
            fieldptr = mxGetField(prhs[1], 0, "ref_signal");
            if(fieldptr){
                custom_ref_sound = true;
                ref_sound_duration_f = mxGetNumberOfElements(fieldptr)/frameSize;
                
                double *ptr = mxGetPr(fieldptr);
                
                for(int i = 0; i<ref_sound_duration_f*frameSize;++i){
                    ref_sound[i]=ptr[i];
                }
                
            }else{
                custom_ref_sound = false;
                
                fieldptr = mxGetField(prhs[1], 0, "ref_freq");
                if(fieldptr){
                    ref_sound_freq = mxGetScalar(fieldptr);
                    if(ref_sound_freq<0) ref_sound_freq=0;
                }else{
                    ref_sound_freq = 200;
                }
                
                fieldptr = mxGetField(prhs[1], 0, "ref_duration");
                if(fieldptr){
                    ref_sound_duration_f = mxGetScalar(fieldptr)*sampleRate/(double)frameSize;
                }else{
                    ref_sound_duration_f = 1*sampleRate/(double)frameSize;
                }
                
                
                fieldptr = mxGetField(prhs[1], 0, "ref_amplitude");
                if(fieldptr){
                    ref_sound_amplitude = mxGetScalar(fieldptr);
                }else{
                    ref_sound_amplitude = 0.2;
                }

            }
            //mexPrintf("p:%i,%i,%i\n",ref_sound_freq,ref_sound_duration_f,ref_sound_amplitude);
            
            fieldptr = mxGetField(prhs[1], 0, "do_var");
            if(fieldptr){
                do_var = bool(mxGetScalar(fieldptr));
            }else{
                do_var = false;
            }
            
            
            fieldptr = mxGetField(prhs[1], 0, "fc");
            if(fieldptr){
                fc = mxGetScalar(fieldptr);
                if(fc<0) fc=0;
            }else{
                fc = 0.01;
            }
           
            fieldptr = mxGetField(prhs[1], 0, "std_dev");
            if(fieldptr){
                std_dev = double(mxGetScalar(fieldptr));
                if(std_dev<0) std_dev=0;
            }else{
                std_dev = 100;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "T_var");
            if(fieldptr){
                T_var_f = int(mxGetScalar(fieldptr));
                if(T_var_f<1) T_var_f=1;
            }else{
                T_var_f = 1;
            }

            fieldptr = mxGetField(prhs[1], 0, "var_quant_size");
            if(fieldptr){
                var_quant_size = int(mxGetScalar(fieldptr));
                if(var_quant_size<0) var_quant_size=0;
            }else{
                var_quant_size = 0;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "do_var_full_trial");
            if(fieldptr){
                do_var_whole_session = bool(mxGetScalar(fieldptr));
            }else{
                do_var_whole_session = false;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "do_control");
            if(fieldptr){
                do_control = bool(mxGetScalar(fieldptr));
            }else{
                do_control = false;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "kp");
            if(fieldptr){
                kp = int(mxGetScalar(fieldptr));
                if(kp<0) kp=0;
            }else{
                kp = 0;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "ki");
            if(fieldptr){
                ki = int(mxGetScalar(fieldptr));
                if(ki<0) ki=0;
            }else{
                ki = 1;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "feedback_gain");
            if(fieldptr){
                feedback_gain = double(mxGetScalar(fieldptr));
                if(feedback_gain<0) feedback_gain=0;
            }else{
                feedback_gain = 1;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "add_pink_noise");
            if(fieldptr){
                add_pink_noise = bool(mxGetScalar(fieldptr));
            }else{
                add_pink_noise = false;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "adaptive_noise_level");
            if(fieldptr){
                adaptive_noise_level = bool(mxGetScalar(fieldptr));
            }else{
                adaptive_noise_level = false;
            }

            fieldptr = mxGetField(prhs[1], 0, "noise_gain");
            if(fieldptr){
                noise_gain = double(mxGetScalar(fieldptr));
                if(noise_gain<0) noise_gain=0;
            }else{
                noise_gain = 0.005;
            }

            fieldptr = mxGetField(prhs[1], 0, "min_noise_level");
            if(fieldptr){
                min_noise_level = double(mxGetScalar(fieldptr))/noise_gain;
                if(min_noise_level<0) min_noise_level=0;
            }else{
                min_noise_level = 0.005/noise_gain;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "max_noise_level");
            if(fieldptr){
                max_noise_level = double(mxGetScalar(fieldptr))/noise_gain;
                if(max_noise_level<min_noise_level) max_noise_level=min_noise_level;
            }else{
                max_noise_level = 0.05/noise_gain;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "start_threshold");
            if(fieldptr){
                start_threshold = double(mxGetScalar(fieldptr));
                if(start_threshold<0) start_threshold=0;
            }else{
                max_noise_level = 0.01;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "stop_threshold");
            if(fieldptr){
                stop_threshold = double(mxGetScalar(fieldptr));
                if(stop_threshold<0) stop_threshold=0;
            }else{
                stop_threshold = 0.01;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "deviceId");
            if(fieldptr){
                deviceId = int(mxGetScalar(fieldptr));
                if(deviceId<0) deviceId=0;
            }else{
                deviceId = 0;
            }
            
            fieldptr = mxGetField(prhs[1], 0, "shifterId");
            if(fieldptr){
                shifterId = int(mxGetScalar(fieldptr));
                if(shifterId<0) shifterId=0;
                if(shifterId>2) shifterId=2;
            }else{
                shifterId = 0;
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
                }
                mxDestroyArray(pink_noise_m);
                mxDestroyArray(num_samples_m);
            }
            
            
            //start stream
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

double sig_power(void){
    if(i_frame < sig_power_width) return 0;
    
    double * ptr = &input_signal[(i_frame+1-sig_power_width)*frameSize];
    
    double sum = 0;
    for(int i=0; i<=sig_power_width*frameSize; ++i){
        sum += ptr[i]*ptr[i];
    }
    
    return sum/(sig_power_width*frameSize);
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
        instrument->noteOn(f, 0.2);

        for(int i = 0; i < ref_sound_duration_f*frameSize; ++i){
            ref_sound[i] = (amp/0.2)*(double) instrument->tick();
            //mexPrintf("%i:%f\n",i,ref_sound[i]);
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
        
        omega = tan(pi*fc/2);
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
    
    if(fc>=1){
        y=x;
    }else{
        y = a0*x+a1*xm1+a2*xm2-b1*ym1-b2*ym2;

        ym2 = ym1;
        ym1 = y;
        xm2 = xm1;
        xm1 = x;
    }
    
    if(var_quant_size > 0){
        y = var_quant_size*round(y/var_quant_size);
    }
    
    return pow(2.0,y/1200.0);
    
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