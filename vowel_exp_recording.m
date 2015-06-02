function data = vowel_exp_recording
    data.subject = 'Linus';
    data.piano_freq = 150;
    save_data = true;
    dir = 'E:\Data\Linus\MATLAB\rec_vowel_exp_data\';
    
    data.condition = 1;
    
    
    switch data.condition
        case 1
            data.condition_name = 'test';
            save_data = false;
            data.mode = 1; %1:shift before voice onset, 2: shift after voice onset
            data.pitch_levels_cents = 200;
            data.num_sessions = 1;
            data.shift_duration_ms = 3*400;
            data.voc_duration_ms = 3000;
            data.shift_onset_interval_ms = 2*[500 500];
            
            data.do_var = 0;
            data.do_var_whole_session = 1;
            data.std_dev = 100;
            data.fc = 0.05*2;

            data.do_control = 0;
            data.kp  = 0;
            data.ki = 0;
            
            data.noise_gain = 0.05;
            
            data.pause_between_sessions_s = 0.7;  
            data.play_ref_whole_session = 1; %1: no, -1: yes
        case 2
            data.condition_name = 'const_shift';
            
            data.mode = 2; %1:shift before voice onset, 2: shift after voice onset
            data.pitch_levels_cents = [-100  0 100];
            data.num_sessions = 60;
            data.shift_duration_ms = 1000;
            data.voc_duration_ms = 2500;
            data.shift_onset_interval_ms = [500 900];
            
            data.do_var = 0;
            data.do_var_whole_session = 0;
            data.std_dev = 0;
            data.fc = 0;

            data.do_control = 0;
            data.kp  = 0;
            data.ki = 0;
            
            data.noise_gain = 0;
            
            data.pause_between_sessions_s = 0.7;  
            data.play_ref_whole_session = 1; %1: no, -1: yes
        case 3
            data.condition_name = 'var_shift';
            
            data.mode = 2; %1:shift before voice onset, 2: shift after voice onset
            data.pitch_levels_cents = [-100  0 100];
            data.num_sessions = 60;
            data.shift_duration_ms = 1000;
            data.voc_duration_ms = 2500;
            data.shift_onset_interval_ms = [500 900];
            
            data.do_var = 1;
            data.do_var_whole_session = 0;
            data.std_dev = 100;
            data.fc = 0.01;

            data.do_control = 0;
            data.kp  = 0;
            data.ki = 0;
            
            data.noise_gain = 0;
            
            data.pause_between_sessions_s = 0.7;  
            data.play_ref_whole_session = 1; %1: no, -1: yes
        case 4
            data.condition_name = 'control_reference';
            
            data.mode = 2; %1:shift before voice onset, 2: shift after voice onset
            data.pitch_levels_cents = 0;
            data.num_sessions = 10;
            data.shift_duration_ms = 0;
            data.voc_duration_ms = 5000;
            data.shift_onset_interval_ms = [4000 4000];
            
            data.do_var = 0;
            data.do_var_whole_session = 0;
            data.std_dev = 0;
            data.fc = 0;

            data.do_control = 0;
            data.kp  = 0;
            data.ki = 0;
            
            data.noise_gain = 0;
            
            data.pause_between_sessions_s = 1.0;  
            data.play_ref_whole_session = 1; %1: no, -1: yes
        case 5
            data.condition_name = 'control';
            
            data.mode = 2; %1:shift before voice onset, 2: shift after voice onset
            data.pitch_levels_cents = 0;
            data.num_sessions = 10;
            data.shift_duration_ms = 0;
            data.voc_duration_ms = 5000;
            data.shift_onset_interval_ms = [4000 4000];
            
            data.do_var = 0;
            data.do_var_whole_session = 0;
            data.std_dev = 0;
            data.fc = 0;

            data.do_control = 1;
            data.kp  = 0;
            data.ki = 2.5;
            
            data.noise_gain = 0;
            
            data.pause_between_sessions_s = 1.0;  
            data.play_ref_whole_session = 1; %1: no, -1: yes
        case 6
            data.condition_name = 'constshort';
            
            data.mode = 2; %1:shift before voice onset, 2: shift after voice onset
            data.pitch_levels_cents = [-100  0 100];
            data.num_sessions = 60;
            data.shift_duration_ms = 400;
            data.voc_duration_ms = 2000;
            data.shift_onset_interval_ms = [500 900];
            
            data.do_var = 0;
            data.do_var_whole_session = 0;
            data.std_dev = 0;
            data.fc = 0;

            data.do_control = 0;
            data.kp  = 0;
            data.ki = 0;
            
            data.noise_gain = 0.002;
            
            data.pause_between_sessions_s = 0.5;  
            data.play_ref_whole_session = 1; %1: no, -1: yes
        case 7
            data.condition_name = 'varshort';
            
            data.mode = 2; %1:shift before voice onset, 2: shift after voice onset
            data.pitch_levels_cents = [-100  0 100];
            data.num_sessions = 60;
            data.shift_duration_ms = 400;
            data.voc_duration_ms = 2000;
            data.shift_onset_interval_ms = [500 900];
            
            data.do_var = 1;
            data.do_var_whole_session = 2;
            data.std_dev = 200;
            data.fc = 0.05;
            
            %double quant_size = 10;
            %int T_var_f = 8;

            data.do_control = 0;
            data.kp  = 0;
            data.ki = 0;
            
            data.noise_gain = 0.002;
            
            data.pause_between_sessions_s = 0.7;  
            data.play_ref_whole_session = 1; %1: no, -1: yes
        otherwise
            fprinf('condition doesn'' exist!\n');
            return;
    end
    
    data.shifter_function = @vowel_shifter_rubberband;
    data.rec_date = datetime('now');
    
    data.frameSize = 64;
    data.Fs = 44100;

    data.pitch_levels = 2.^(data.pitch_levels_cents/1200);
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
        data.shifter_function(data.mode, data.pitch_level_sqs(i), data.voc_duration_f, data.play_ref_whole_session*data.piano_freq, data.shift_onset_f(i), data.shift_duration_f, data.do_var, data.std_dev, data.fc, data.do_control, data.kp, data.ki, data.do_var_whole_session, data.noise_gain);
        while(data.shifter_function(0) == 0)
            pause(0.2);
        end
        [data.y_r{i}, data.y_ps{i}, data.voice_onset_f(i),data.static_pitch_factor_sqs{i}, data.var_pitch_factor_sqs{i}, data.control_pitch_factor_sqs{i},data.detected_pitch{i}] = data.shifter_function(-1);
        pause(data.pause_between_sessions_s);
    end
    data.voice_onset_ms = data.voice_onset_f * data.frameSize * 1000 / data.Fs;
    
    if save_data
        [~,M,D] = datevec(data.rec_date);
        filename = sprintf('%s%s%02i%02i_%s.mat',dir,data.subject,M,D,data.condition_name);
        [~,name,~] = fileparts(filename);
        if exist(filename,'file')
            choice = questdlg(sprintf('The file %s already exists. Do you want to overwrite it?',name),'Warning!', 'Yes', 'No', 'No');
        else
            choice = 'Yes';
        end
        if strcmpi(choice,'Yes')
            fprintf('save file: %s\n',name);
            save(filename,'data');
        end
    end

end