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
const double minlevelduration_s = 1;
const int minlevelduration_f = (int) round(minlevelduration_s*sampleRate/(double)frameSize);

const double minlowampduration_s = 0.01;
const int minlowampduration_f = (int) ceil(minlowampduration_s*sampleRate/(double)frameSize);

const float threshold = 0.005;
const float std_dev = 10;
const float delta = 0.01;
const float p_is_var = 1;

double pitch_levels[5] = {-80.0,-40.0,40.0,80.0};
                           
default_random_engine generator;
uniform_int_distribution<int> distribution_lvl_choice(0,3);
bernoulli_distribution distribution_is_var(p_is_var);
normal_distribution<double> distribution_pitch_var(0.0,std_dev);

float *output;
int i_frame;

int i_session;


int i_low_amp;

const size_t data_array_length = 360*44100;
float input_signal[data_array_length];
float output_signal[data_array_length];
int levels[1000];
int session_starts[1000];

RubberBandStretcher::Options options =  RubberBandStretcher::OptionProcessRealTime |
                                        RubberBandStretcher::OptionFormantPreserved |
                                        RubberBandStretcher::OptionWindowShort |
                                        RubberBandStretcher::OptionPitchHighConsistency;
RubberBandStretcher shifter(sampleRate,1,options,1.0,pow(2.0,0/12.0));

//variables only used in the tick function
bool is_var;
double var;

void init_global_vars(void){
    i_frame = -1;
    i_session = 0;
    i_low_amp = 0;
    
    shifter.reset();
    shifter.setPitchScale(1.0);
            
    session_starts[0] = 0;
    levels[0] = -1;
    
    is_var = false;
    var = 0;
    
    memset(input_signal,0,sizeof(input_signal));
    memset(output_signal,0,sizeof(output_signal));
    memset(levels,0,sizeof(levels));
    memset(session_starts,0,sizeof(session_starts));
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
    
    
    if(i_low_amp >= minlowampduration_f && i_frame - session_starts[i_session] >= minlevelduration_f) {
        i_session++;
        
        is_var = distribution_is_var(generator);
        var = 0;
        
        session_starts[i_session] = i_frame;
        
        if(i_session % 2 == 1){
            levels[i_session] = distribution_lvl_choice(generator);
            shifter.setPitchScale(pow(2.0,(pitch_levels[levels[i_session]]+var)/1200.0));
        }else{
            levels[i_session] = -1;
            shifter.setPitchScale(1.0);
        }
        
        
        
        
        //int lat_sw = shifter.getLatency();
        //mexPrintf("level: %i isvar: %i lat: %i\n", lvl, is_var, lat_sw);
        //mexPrintf("amp: %f\n", amp);
        mexPrintf("level: %i session: %i\n", levels[i_session], i_session);
    }
    
//         if(is_var){
//             var = (1-delta) * var + delta * distribution_pitch_var(generator);
//         }
    
    //var= 0;
    //lvl = 2;
    //shifter.setPitchScale(pow(2.0,(pitch_levels[lvl]+var)/1200.0));
    //shifter.setPitchScale(1.0);
    
    shifter.process(&ibuffer, nBufferFrames,0);
    
    //int avl = shifter.RubberBandStretcher::available();
    
//     if(avl>nBufferFrames){
//         //int u = shifter.retrieve(&output, avl-nBufferFrames);
//         //mexPrintf("a: %i frame: %i var: %i session: %i\n", u, i_frame,lvl,i_session);
//     }
//     else if(avl<nBufferFrames){
//         //mexPrintf("r: %i frame: %i var: %i session: %i\n", avl, i_frame,lvl,i_session);
//     }
    
    int r = shifter.retrieve(&output, nBufferFrames);
    
    memcpy(&output_signal[i_frame*frameSize],output,r*sizeof(float));
    
//     if(r!=nBufferFrames){
//          mexPrintf("r: %i frame: %i lvl: %i session: %i\n", r, i_frame,lvl,i_session);
//     }

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

	RtAudio dac;
    
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
    mexCallMATLAB(0, NULL, 0, NULL, "pause");
    
    try {
		dac.closeStream();
	}
	catch ( RtAudioError &error ) {
		error.printMessage();
        return;
	}
    
    delete[] output;
    
    return;
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
    if(nlhs>4) {
        mexErrMsgIdAndTxt("rt_pitch_shifter:nlhs", "only 4 outputs allowed.");
    }
    
    init_global_vars();
    rt_pitch_shifter();
    
    mxArray *input_signal_m = mxCreateDoubleMatrix((i_frame+1)*frameSize,1,mxREAL);
    double *input_signal_mp = mxGetPr(input_signal_m);
    mxArray *output_signal_m = mxCreateDoubleMatrix((i_frame+1)*frameSize,1,mxREAL);
    double *output_signal_mp = mxGetPr(output_signal_m);
    
    for(int i = 0; i<(i_frame+1)*frameSize; ++i){
        input_signal_mp[i] = (double) input_signal[i];
        output_signal_mp[i] = (double) output_signal[i];
    }
    
    mxArray *session_starts_m = mxCreateDoubleMatrix(i_session+1,1,mxREAL);
    double *session_starts_mp = mxGetPr(session_starts_m);
    mxArray *levels_m = mxCreateDoubleMatrix(i_session+1,1,mxREAL);
    double *levels_mp = mxGetPr(levels_m);
    
    for(int i = 0; i<i_session+1; ++i){
        session_starts_mp[i] = (double) session_starts[i]*frameSize+1;
        levels_mp[i] = (double) levels[i];
    }

    if(nlhs >= 1) plhs[0] = input_signal_m;
    if(nlhs >= 2) plhs[1] = output_signal_m;
    if(nlhs >= 3) plhs[2] = session_starts_m;
    if(nlhs >= 4) plhs[3] = levels_m;
    
    return;
}

float frameAmplitude(float *samples) {
    float sum = 0;
    for(int i = 0; i < frameSize; ++i) {
        sum += abs(samples[0]);
    }
    return sum/frameSize;
        
}