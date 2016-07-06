function data = perc_adapt_exp_recording()
% pitch shift experiment according to Lindner et al. 2008 

    data.subject = 'Richard';
    
    data.ref_freq = 200;
    data.noise_gain = 0;
    
    data.voc_duration = 2;
    
    data.random_start_shift = true; %random start shift for endless loop playback
    
    data.condition = 1;
    
    % final pitch shift
    data.pitch_shift_cents = -300;
    
    % pitch shift maybe noisy
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
    
    % reference pitch maybe noisy too
    data.ref_shifts_cents = [-150 -50 50 150];
    data.ref_freqs = 2.^(data.ref_shifts_cents/1200) * data.ref_freq;
    
    % there are four kind of trials and three phases (before shift trials, transition
    % trials, after shift trials).
    
    % Define the number of trials in each phase here
    data.n_pre_fb_trials = 1;
    data.n_pre_mct_trials = 0;
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
    
    if data.random_start_shift
        data.startShift = (2*rand(1,data.n_trials)-1)*600;
    end
    
    % get the random sequence of pitch shifts and reference frequencies
    for i=1:data.n_trials
        data.pitch_level_var_sqs_cents(i) = data.pitch_level_sqs_cents(i) + data.var_mult * data.var_pitch_shift_cents(ceil(length(data.var_pitch_shift_cents)*rand(1)));
        data.ref_freq_sqs(i) = data.ref_freqs(ceil(length(data.ref_freqs)*rand(1)));
    end
    
    % get order of trials
    
    % trial order in phase 1 and 3 is random
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
    
    
    %% setup gui
    gui = figure;
    gui.MenuBar = 'none';
    
    %status text
    txt_sn = uicontrol(gui, 'Style','text',...
            'Units', 'normalized',...
            'Position',[0 0.7 1 0.25],...
            'String','',...
            'FontSize',14);
        
    %text to signal trial end 
    txt_stop = uicontrol(gui, 'Style','text',...
            'Units', 'normalized',...
            'Position',[0 0.5 1 0.25],...
            'String','Stop!',...
            'ForegroundColor', 'r',...
            'Visible', 'off',...
            'FontSize',14);
    
    %% prepare PsychPitchShifter settings
    params.shifterId = 2;
    params.windowSize = 2;
    params.deviceId = 0;
    params.shift_full_trial = true;
    params.voc_duration = data.voc_duration;
    params.start_threshold = 0.01;
    params.stop_threshold = 0.005;
    params.shift_full_trial = 1;
    
    params.volume_normalization = 1;
    params.tau = 0.3;
    params.delta = 0.01;

    %main loop    
    for i=1:data.n_trials
        txt_stop.Visible = 'off';
        
        %trial dependent settings
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
        
        %start recording
        PsychPitchShifter(1, params);
        drawnow;
        
        %wait till end of trial
        while(PsychPitchShifter(0) == 0)
            pause(0.2);
        end
        
        %wait till voice offset
        if PsychPitchShifter(0) == 1
            txt_stop.Visible = 'on'; %signal trial end
            while(PsychPitchShifter(0) ~= 2)
                pause(0.2);
            end
        end
        % stop recording
        [data.y_r{i}, data.y_ps{i}, data.voice_onset_s(i),~, ~, ~,data.detected_pitch{i}] = PsychPitchShifter(-1);
       
        d = data.y_r{i}(round((data.voice_onset_s(i)+0.2)*data.Fs):round((data.voice_onset_s(i)+data.voc_duration-0.2)*data.Fs));
       
        % let subject guess shift
        if data.random_start_shift
            [data.perceived_produced_pitch(i),data.invalid_trials(i)] = get_perceived_shift(d,data.startShift(i),params);
        else
            [data.perceived_produced_pitch(i),data.invalid_trials(i)] = get_perceived_shift(d,data.pitch_level_var_sqs_cents(i),params);
        end
    end
    
    delete(gui);

end