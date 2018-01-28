params.shifterId = 2;
params.pitch_factor = 2^(50/1200);
params.voc_duration = 2;
params.start_threshold = 0.01;
params.shift_duration = 1;
params.windowSize = 1;
params.rubberbandOptionPitch = 2;
params.rubberbandOptionPhase = 1;

fs = 44100;
frameSize = 64;


T=4;
T_f = floor(T*fs/frameSize);
T= T_f*frameSize/fs;

t=(1/fs:1/fs:T);

f=600;
sig = sin(2*pi*f*t)+sin(2*pi*2*f*t)+sin(2*pi*3*f*t);

load tests/gendata/flat_female_c1300_HPF.mat

%%
PsychPitchShifter(2,params,realHPF);
[x,y,voice_on,stc_pf,var_pf,ctrl_pf,dpitch] = PsychPitchShifter(-1);

soundsc(y,fs);