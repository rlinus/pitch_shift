function data = discrete_exp_recording(data)
    data.subject = 'Richard';
    
    data.guess_shift = 1;
    
    data.piano_freq = 0;
    
    data.voc_duration_ms = 2000;
    data.std = 0;
    
    data.noise_gain = 0.000;
    data.pause_between_sessions_s = 0.7;
    
% %%
%     
%     data.pitch_levels_cents = 3*[-100 -50 0 50 100];
%     data.num_sessions = 10;
%     
%     data.pitch_level_sqs_cents = zeros(1,data.num_sessions);
%     for i=1:data.num_sessions
% 
%         lvl = ceil(length(data.pitch_levels_cents)*rand(1));
%         data.pitch_level_sqs_cents(i) = data.pitch_levels_cents(lvl);
% 
%     end
%     data.pitch_level_var_sqs_cents = data.pitch_level_sqs_cents + data.std*randn(1,data.num_sessions);
%     
%     data.pitch_level_var_sqs = 2.^(data.pitch_level_var_sqs_cents/1200);
%     
%     %data.startShift = (2*rand(1,data.num_sessions)-1)*(max(data.pitch_levels_cents)-min(data.pitch_levels_cents))/2+(max(data.pitch_levels_cents)+min(data.pitch_levels_cents))/2;
%     %data.startShift = (2*rand(1,data.num_sessions)-1)*400;
    
    
%    
    
    data.n_pre_trials = 12;
    data.n_post_trials = 8;
    data.n_trans_trials = 12;
    data.num_sessions = data.n_pre_trials + data.n_post_trials + data.n_trans_trials;
    
     
    
    data.pitch_shift_cents = 300*sign(randn(1));
    
    
    data.pitch_levels_cents = [0 data.pitch_shift_cents];

    inc = data.pitch_shift_cents/(data.n_trans_trials+1);
    data.pitch_level_sqs_cents = [zeros(1,data.n_pre_trials),(inc:inc:data.pitch_shift_cents-inc), ones(1,data.n_post_trials)*data.pitch_shift_cents];
    
    data.pitch_level_var_sqs_cents = data.pitch_level_sqs_cents + data.std*randn(1,data.num_sessions);
    
    data.pitch_level_var_sqs = 2.^(data.pitch_level_var_sqs_cents/1200);

%%
    

    data.frameSize = 64;
    data.Fs = 44100;
    
    data.shifter_function = @vowel_shifter_rubberband;
    data.rec_date = datetime('now');
    
    data.voc_duration_f = round(data.voc_duration_ms * data.Fs / data.frameSize / 1000);
    data.voc_duration_ms = data.voc_duration_f * data.frameSize * 1000 / data.Fs;
    
    for i=1:data.num_sessions
        fprintf('session %i...\n',i);
%         %reset trial
%          data.shifter_function(1, 1, data.voc_duration_f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, data.noise_gain);
%         while(data.shifter_function(0) == 0)
%             pause(0.2);
%         end
%         data.shifter_function(-1);
        
        
        %probe trial
        data.shifter_function(1, data.pitch_level_var_sqs(i), data.voc_duration_f, data.piano_freq, 0, 0, 0, 0, 0, 0, 0, 0, 0, data.noise_gain);
        while(data.shifter_function(0) == 0)
            pause(0.2);
        end
        [data.y_r{i}, data.y_ps{i}, data.voice_onset_f(i),data.static_pitch_factor_sqs{i}, data.var_pitch_factor_sqs{i}, data.control_pitch_factor_sqs{i},data.detected_pitch{i}] = data.shifter_function(-1);
        data.voice_onset_ms(i) = data.voice_onset_f(i) * data.frameSize * 1000 / data.Fs;
        if data.guess_shift
            d = data.y_r{i}(data.voice_onset_f(i) * data.frameSize+0.2*data.Fs:(data.voice_onset_f(i)+data.voc_duration_f) * data.frameSize-0.1*data.Fs);
            %data.perceived_produced_pitch(i) = get_perceived_shift(d,data.startShift(i));
            data.perceived_produced_pitch(i) = get_perceived_shift(d,data.pitch_level_var_sqs_cents(i));
        else
            pause(data.pause_between_sessions_s);
        end
    end
    

end

