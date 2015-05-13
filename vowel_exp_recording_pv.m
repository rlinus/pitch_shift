function data = vowel_exp_recording_pv  
    data.frameSize = 64;
    data.Fs = 44100;
    
    data.mode = 2; %1:shift before voice onset, 2: shift after voice onset

    data.pitch_levels_cents = [-100 0 100];%[-60 -30 0 30 60];
    data.pitch_levels = 2.^(data.pitch_levels_cents/1200);

    data.num_sessions = 5;
    data.shift_duration_ms = 800;
    data.voc_duration_ms = 2000;
    data.shift_onset_interval_ms = [400 800];
    
    data.piano_freq = 200;
    data.play_ref_whole_session = 1; %1: no, -1: yes
    
    data.do_var = 0;
    data.std_dev = 100;
    data.fc = 0.005;
    
    data.do_control = 0;
    data.kp  = 0;
    data.ki = 5;
    
    data.pause_between_sessions_s = 0.5;
    
    
    
    data.rec_date = datetime('now');

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
        vowel_shifter_pv(data.mode, data.pitch_level_sqs(i), data.voc_duration_f, data.play_ref_whole_session*data.piano_freq, data.shift_onset_f(i), data.shift_duration_f, data.do_var, data.std_dev, data.fc, data.do_control, data.kp, data.ki);
        while(vowel_shifter_pv(0) == 0)
            pause(0.2);
        end
        [data.y_r{i}, data.y_ps{i}, data.voice_onset_f(i),data.static_pitch_factor_sqs{i}, data.var_pitch_factor_sqs{i}, data.control_pitch_factor_sqs{i},data.detected_pitch{i}] = vowel_shifter_pv(-1);
        pause(data.pause_between_sessions_s);
    end
    data.voice_onset_ms = data.voice_onset_f * data.frameSize * 1000 / data.Fs;

end