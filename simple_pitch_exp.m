Fs = 44100;

f0_low=95; % take only pitch values above that level
f0_high=300;%140; % below that level
max_pitch_change=5; % disregard pitch changes (from one point to the next) by more than 10.

rng('shuffle');
pitch_cents = round((rand(1,1)-0.5)*400);
pitch_factor = 2^(pitch_cents/1200);

num_intro_frames = 2*2100;

num_add_frames = round(rand(1,1)*1400);

num_frames = num_intro_frames+num_add_frames;

wait_time = num_frames*64/44100+6;

rt_pitch_shifter(num_frames,pitch_factor);
pause(wait_time);
[is,os] = rt_pitch_shifter(-1);


pitch_change_time_ms = num_frames*64/44.1;
[f0_time_is,f0_value_is,SHR_is,f0_candidate_is]=shrp(is,Fs);
[f0_time_os,f0_value_os,SHR_os,f0_candidate_os]=shrp(os,Fs);

cand=abs(diff(f0_value_is))<max_pitch_change & f0_value_is(2:end)>f0_low & f0_value_is(2:end)<f0_high;
cand = [logical(0); cand];
f0nan_is=nan(size(f0_value_is)); f0nan_is(cand)=f0_value_is(cand);

mean_before = mean(f0_value_is(f0_time_is<pitch_change_time_ms & cand));
mean_afterwards = mean(f0_value_is(f0_time_is>pitch_change_time_ms+500 & cand));
pitch_change = 1200 * log2(mean_afterwards/mean_before);

clf;
plot(f0_time_is,f0_value_is);
hold on;
plot(f0_time_is,f0nan_is,'r');
%plot(f0_time_os,f0_value_os);
plot([pitch_change_time_ms, pitch_change_time_ms],[0,max(f0_value_is)]);

