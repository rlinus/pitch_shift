function data = discrete_exp_recording(data)
    data.subject = 'Linus';
    data.piano_freq = 150;
    
    data.voc_duration_ms = 2000;
    data.voc_duration_f = round(data.voc_duration_ms * data.Fs / data.frameSize / 1000);
    data.voc_duration_ms = data.voc_duration_f * data.frameSize * 1000 / data.Fs;
    
    data.n_pre_trials = 2*10;
    data.n_post_trials = 2*10;
    data.n_trans_trials = 4;
    data.num_sessions = data.n_pre_trials + data.n_post_trials + data.n_trans_trials;
    
    data.noise_gain = 0.002;
    data.pause_between_sessions_s = 0.7; 
    
    data.pitch_shift_cents = 100*sign(randn(1));
    data.std = 25;
    
    data.pitch_levels_cents = [0 data.pitch_shift_cents];
    
    data.shifter_function = @vowel_shifter_rubberband;
    data.rec_date = datetime('now');
    
    data.frameSize = 64;
    data.Fs = 44100;
    
    inc = data.pitch_shift_cents/(data.n_trans_trials+1);
    data.pitch_level_sqs_cents = [zeros(1,data.n_pre_trials),(inc:inc:data.pitch_shift_cents-inc), ones(1,data.n_post_trials)*data.pitch_shift_cents];
    
    data.pitch_level_var_sqs_cents = data.pitch_level_sqs_cents + data.std*randn(1,data.num_sessions);
    
    data.pitch_level_var_sqs = 2.^(data.pitch_level_var_sqs_cents/1200);
    
    for i=1:data.num_sessions
        fprintf('session %i...\n',i);
        data.shifter_function(1, data.pitch_level_var_sqs(i), data.voc_duration_f, data.play_ref_whole_session*data.piano_freq, 0, 0, 0, 0, 0, 0, 0, 0, 0, data.noise_gain);
        while(data.shifter_function(0) == 0)
            pause(0.2);
        end
        [data.y_r{i}, data.y_ps{i}, data.voice_onset_f(i),data.static_pitch_factor_sqs{i}, data.var_pitch_factor_sqs{i}, data.control_pitch_factor_sqs{i},data.detected_pitch{i}] = data.shifter_function(-1);
        pause(data.pause_between_sessions_s);
    end
    data.voice_onset_ms = data.voice_onset_f * data.frameSize * 1000 / data.Fs;

end

