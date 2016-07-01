function data = perc_adapt_exp_recording()
    data.subject = 'Richard';
    
    data.ref_freq = 200;
    data.noise_gain = 0;
    
    data.voc_duration_ms = 2*2000;
    
    data.random_start_shift = true;
    
    data.condition = 1;
    
    data.pitch_shift_cents = -300;%*sign(randn(1));
    
    data.var_pitch_shift_cents = [-100 -50 -25 0 25 50 100];
    
    switch data.condition
        case 1
            data.condition_name = 'no_var';
            data.var_mult = 0;
        case 2
            data.condition_name = 'small_discrete_var';
            data.var_mult = 1;
        case 3
            data.condition_name = 'big_discrete_var';
            data.var_mult = 2;
        otherwise
            warning('Condition not defined!');
            return;
    end
    
    data.ref_shifts_cents = [-150 -50 50 150];
    data.ref_freqs = 2.^(data.ref_shifts_cents/1200) * data.ref_freq;
    
    data.n_pre_fb_trials = 0;
    data.n_pre_mct_trials = 0;
    data.n_pre_pct_trials = 1;
    
    data.n_post_fb_trials = 0; 
    data.n_post_mct_trials = 0;
    data.n_post_pct_trials = 0;
    
    data.n_trans_trials = 0;
    
    data.n_pre_trials = data.n_pre_fb_trials + data.n_pre_mct_trials + data.n_pre_pct_trials;
    data.n_post_trials = data.n_post_fb_trials + data.n_post_mct_trials + data.n_post_pct_trials;
    data.n_trials = data.n_pre_trials + data.n_post_trials + data.n_trans_trials;
    
    inc = data.pitch_shift_cents/(data.n_trans_trials+1);
    data.pitch_level_sqs_cents = [zeros(1,data.n_pre_trials),(inc:inc:data.pitch_shift_cents-inc), ones(1,data.n_post_trials)*data.pitch_shift_cents];
    
    if data.random_start_shift
        data.startShift = (2*rand(1,data.n_trials)-1)*600;
    end
    
    for i=1:data.n_trials
        data.pitch_level_var_sqs_cents(i) = data.pitch_level_sqs_cents(i) + data.var_mult * data.var_pitch_shift_cents(ceil(length(data.var_pitch_shift_cents)*rand(1)));
        data.ref_freq_sqs(i) = data.ref_freqs(ceil(length(data.ref_freqs)*rand(1)));
    end
    
    pre_kind = [1*ones(1,data.n_pre_fb_trials), 2*ones(1,data.n_pre_mct_trials), 3*ones(1,data.n_pre_pct_trials)];
    pre_kind = pre_kind(randperm(data.n_pre_trials));
    post_kind = [1*ones(1,data.n_post_fb_trials), 2*ones(1,data.n_post_mct_trials), 3*ones(1,data.n_post_pct_trials)];
    post_kind = post_kind(randperm(data.n_post_trials));
    
    data.kind_of_trials = [pre_kind 1*ones(1,data.n_trans_trials) post_kind];
    
    data.pitch_level_var_sqs = 2.^(data.pitch_level_var_sqs_cents/1200);
    
    data.invalid_trials = zeros(data.n_trials,1);
    
    data.frameSize = 64;
    data.Fs = 44100;

    data.rec_date = datetime('now');
    
    data.voc_duration_f = round(data.voc_duration_ms * data.Fs / data.frameSize / 1000);
    data.voc_duration_ms = data.voc_duration_f * data.frameSize * 1000 / data.Fs;
    
    gui = figure;
    gui.MenuBar = 'none';
    
    txt_sn = uicontrol(gui, 'Style','text',...
            'Units', 'normalized',...
            'Position',[0 0.7 1 0.25],...
            'String','',...
            'FontSize',14);
        
    txt_stop = uicontrol(gui, 'Style','text',...
            'Units', 'normalized',...
            'Position',[0 0.5 1 0.25],...
            'String','Stop!',...
            'ForegroundColor', 'r',...
            'Visible', 'off',...
            'FontSize',14);
        
    params.shifterId = 0;
    params.deviceId = 0;
    params.shift_full_trial = true;
    params.voc_duration = data.voc_duration_ms/1000;
    params.shift_duration = data.shift_duration_ms/1000;
    params.start_threshold = 0.1;
    params.stop_threshold = 0.1;

        
    for i=1:data.n_trials
        txt_stop.Visible = 'off';
        
        params.pitch_factor = data.pitch_level_var_sqs(i);
        
        switch data.kind_of_trials(i)
            case 1 %feedback trial
                params.add_pink_noise = 0;
                params.feedback_gain = 1;
                params.play_ref_sound = 0;

                txt_sn.String = sprintf('Session %i: feedback trial', i);
            case 2 %motor catch trial
                params.add_pink_noise = 1;
                params.feedback_gain = 0;
                params.play_ref_sound = 1;
                params.ref_freq = data.ref_freq_sqs(i);
                params.ref_amplitude = 0.1;
                
                txt_sn.String = sprintf('Session %i: motor catch trial', i);
            case 3 %perceptual catch trial
                params.add_pink_noise = 1;
                params.feedback_gain = 0;
                params.play_ref_sound = 0;
                
                txt_sn.String = sprintf('Session %i: perceptual catch trial', i);
        end
        PsychPitchShifter(1, params);
        drawnow;
        while(pitch_shifter(0) == 0)
            pause(0.2);
        end
        if PsychPitchShifter(0) == 1
            txt_stop.Visible = 'on';
            while(PsychPitchShifter ~= 2)
                pause(0.2);
            end
        end
        [data.y_r{i}, data.y_ps{i}, data.voice_onset_f(i),~, ~, ~,data.detected_pitch{i}] = PsychPitchShifter(-1);
        data.voice_onset_ms(i) = data.voice_onset_f(i) * data.frameSize * 1000 / data.Fs;
       
        d = data.y_r{i}(data.voice_onset_f(i) * data.frameSize+0.2*data.Fs:(data.voice_onset_f(i)+data.voc_duration_f) * data.frameSize-0.1*data.Fs);
       
        if data.random_start_shift
            [data.perceived_produced_pitch(i),data.invalid_trials(i)] = get_perceived_shift(d,data.startShift(i));
        else
            [data.perceived_produced_pitch(i),data.invalid_trials(i)] = get_perceived_shift(d,data.pitch_level_var_sqs_cents(i));
        end
    end
    
    delete(gui);

end