function data = vowel_exp_recording_pv  
    data.frameSize = 64;
    data.Fs = 44100;
    
    data.mode = 2; %1:shift before voice onset, 2: shift after voice onset

    data.pitch_levels_cents = [10];%[-60 -30 0 30 60];
    data.pitch_levels = 2.^(data.pitch_levels_cents/1200);

    data.num_sessions = 1;
    data.shift_duration_ms = 2000;
    data.voc_duration_ms = 4000;
    
    data.piano_freq = 300;
    data.play_ref_whole_session = 1; %1: no, -1: yes
    
    data.rec_date = datetime('now');

    data.shift_onset_interval_ms = [1000 1000]; %[400 800];
    data.shift_onset_ms = min(data.shift_onset_interval_ms) + rand(data.num_sessions,1) * abs(diff(data.shift_onset_interval_ms));
    
    data.pitch_level_sqs_cents = zeros(data.num_sessions, 1);
    for i=1:data.num_sessions

        lvl = ceil(length(data.pitch_levels_cents)*rand(1));
        data.pitch_level_sqs_cents(i) = data.pitch_levels_cents(lvl);

    end
    data.pitch_level_sqs = 2.^(data.pitch_level_sqs_cents/1200);

    data.shift_onset_f = round(data.shift_onset_ms * data.Fs / data.frameSize / 1000);
    data.shift_duration_f = round(data.shift_duration_ms * data.Fs / data.frameSize / 1000);
    data.voc_duration_f = round(data.voc_duration_ms * data.Fs / data.frameSize / 1000);
    
    data.shift_onset_ms = data.shift_onset_f * data.frameSize * 1000 / data.Fs;
    data.shift_duration_ms = data.shift_duration_f * data.frameSize * 1000 / data.Fs;
    data.voc_duration_ms = data.voc_duration_f * data.frameSize * 1000 / data.Fs;

    for i=1:data.num_sessions
        fprintf('session %i...\n',i);
        vowel_shifter_pv(data.mode, data.pitch_level_sqs(i), data.voc_duration_f, data.play_ref_whole_session*data.piano_freq, data.shift_onset_f(i), data.shift_duration_f);
        while(vowel_shifter_pv(0) == 0)
            pause(0.2);
        end
        [data.y_r{i}, data.y_ps{i}, data.voice_onset_f(i),data.samples_available{i},data.sw_latency{i}] = vowel_shifter_pv(-1);
        pause(0.3);
    end
    data.voice_onset_ms = data.voice_onset_f * data.frameSize * 1000 / data.Fs;
    
    data.sessions_with_holes = [];
    for i=1:data.num_sessions
        holes = find(data.samples_available{i}<data.frameSize);
        holes = holes(holes>=data.voice_onset_f(i));
        if ~isempty(holes)
            fprintf('session %i has holes after voice onset (first hole at frame %i).\n)',i,holes(1));
            data.holes{i} = holes;
            data.sessions_with_holes(end+1) = i;
        end
    end
end


% data.freq_range = [50 300];
% [f0_time,f0_value]=shrp(y_r,Fs,freq_range);
% [f0_time2,f0_value2]=shrp(y_ps,Fs,freq_range);
% 
% plot(f0_time,f0_value,'b');
% hold on
% plot(f0_time2,f0_value2,'g');
% plot([voice_onset_ms voice_onset_ms], freq_range, 'r');
% plot([voice_onset_ms+shift_onset_ms voice_onset_ms+shift_onset_ms], freq_range, 'r');