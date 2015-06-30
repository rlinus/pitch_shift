function data = perc_adapt_exp_recording()
    data.subject = 'Gagan';
    
    data.ref_freq = 200;
    data.noise_gain = 0.10;
    
    data.voc_duration_ms = 2000;
    
    data.condition = 1;
    
    data.pitch_shift_cents = 300*sign(randn(1));
    
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
    data.n_pre_mct_trials = 5;
    data.n_pre_pct_trials = 0;
    
    data.n_post_fb_trials = 0;
    data.n_post_mct_trials = 0;
    data.n_post_pct_trials = 0;
    
    data.n_trans_trials = 0;
    
    data.n_pre_trials = data.n_pre_fb_trials + data.n_pre_mct_trials + data.n_pre_pct_trials;
    data.n_post_trials = data.n_post_fb_trials + data.n_post_mct_trials + data.n_post_pct_trials;
    data.n_trials = data.n_pre_trials + data.n_post_trials + data.n_trans_trials;
    
    inc = data.pitch_shift_cents/(data.n_trans_trials+1);
    data.pitch_level_sqs_cents = [zeros(1,data.n_pre_trials),(inc:inc:data.pitch_shift_cents-inc), ones(1,data.n_post_trials)*data.pitch_shift_cents];
    
    data.startShift = (2*rand(1,data.n_trials)-1)*600;
    
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
    
    data.frameSize = 64;
    data.Fs = 44100;
    
    data.shifter_function = @vowel_shifter_rubberband;
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
        
    for i=1:data.n_trials
        txt_stop.Visible = 'off';
        
        switch data.kind_of_trials(i)
            case 1 %feedback trial
            	data.shifter_function(1, data.pitch_level_var_sqs(i), data.voc_duration_f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1);
                txt_sn.String = sprintf('Session %i: feedback trial', i);
            case 2 %motor catch trial
            	data.shifter_function(1, data.pitch_level_var_sqs(i), data.voc_duration_f, data.ref_freq_sqs(i), 0, 0, 0, 0, 0, 0, 0, 0, 0, data.noise_gain,0);
                txt_sn.String = sprintf('Session %i: motor catch trial', i);
            case 3 %perceptual catch trial
                data.shifter_function(1, data.pitch_level_var_sqs(i), data.voc_duration_f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, data.noise_gain,0);
                txt_sn.String = sprintf('Session %i: perceptual catch trial', i);
        end
        drawnow;
        while(data.shifter_function(0) == 0)
            pause(0.2);
        end
        if data.shifter_function(0) == 1
            txt_stop.Visible = 'on';
            while(data.shifter_function(0) ~= 2)
                pause(0.2);
            end
        end
       [data.y_r{i}, data.y_ps{i}, data.voice_onset_f(i),~, ~, ~,data.detected_pitch{i}] = data.shifter_function(-1);
       data.voice_onset_ms(i) = data.voice_onset_f(i) * data.frameSize * 1000 / data.Fs;
       
       d = data.y_r{i}(data.voice_onset_f(i) * data.frameSize+0.2*data.Fs:(data.voice_onset_f(i)+data.voc_duration_f) * data.frameSize-0.1*data.Fs);
       data.perceived_produced_pitch(i) = get_perceived_shift(d,data.startShift(i));
       %data.perceived_produced_pitch(i) = get_perceived_shift(d,data.pitch_level_var_sqs_cents(i));
    end
    
    delete(gui);

end