params.shifterId = 1;
params.pitch_factor = 0.9;
params.voc_duration = 2;

sampleRate = 44100;
frameSize = 64;

pitchShift=2^(0/1200);

T=4;
T_f = floor(T*sampleRate/frameSize);
T= T_f*frameSize/sampleRate;

t=(1/sampleRate:1/sampleRate:T);

f=600;
sig = sin(2*pi*f*t);

PsychPitchShifter(2,params,sig);
[x,y,voice_on,stc_pf,var_pf,ctrl_pf,dpitch] = PsychPitchShifter(-1);