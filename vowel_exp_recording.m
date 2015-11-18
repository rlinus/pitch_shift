function [data,params] = vowel_exp_recording
    data.subject = 'Linus2';
    data.piano_freq = 200;
    save_data = true;
    dir = 'C:\Users\Linus\Documents\MATLAB\test_rec\';
    
    data.condition = 1;
    
    
    switch data.condition
        case 1
            data.condition_name = 'test';
            save_data = false;
            data.mode = 2; %1:shift before voice onset, 2: shift after voice onset
            data.pitch_levels_cents = 200;
            data.num_sessions = 1;
            data.shift_duration_ms = 1000;
            data.voc_duration_ms = 3000;
            data.shift_onset_interval_ms = 2*[500 500];
            
            data.do_var = 0;
            data.do_var_whole_session = 0;
            data.std_dev = 100;
            data.fc = 0.05*2;

            data.do_control = 0;
            data.kp  = 0;
            data.ki = 0;
            
            data.noise_gain = 0.0;
            
            data.pause_between_sessions_s = 0.7;  
            data.play_ref_whole_session = -1; %1: no, -1: yes
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
    
    params.shifterId = 0;
    params.deviceId = 0;
    params.shift_after_voice_onset = data.mode==2;
    params.voc_duration = data.voc_duration_ms/1000;
    params.shift_duration = data.shift_duration_ms/1000;
    params.play_ref_sound = 1;
    params.ref_freq = data.piano_freq;
    params.ref_amplitude = 0.1;
    params.add_pink_noise = data.noise_gain ~=0;
    params.noise_gain = data.noise_gain;
    params.do_var = data.do_var;
    params.do_var_whole_session = data.do_var_whole_session;
    params.std_dev = data.std_dev;
    params.fc = data.fc;
    params.do_control = data.do_control;
    params.ki = data.ki;
    params.kp = data.kp;
    params.feedback_gain = 1;
    

    for i=1:data.num_sessions
        fprintf('session %i...\n',i);
        params.shift_onset = data.shift_onset_ms(i)/1000;
        params.pitch_factor = data.pitch_level_sqs(i);
        pitch_shifter(1, params);
        while(pitch_shifter(0) == 0)
            pause(0.2);
        end
        [data.y_r{i}, data.y_ps{i}, data.voice_onset_s(i),data.static_pitch_factor_sqs{i}, data.var_pitch_factor_sqs{i}, data.control_pitch_factor_sqs{i},data.detected_pitch{i}] = pitch_shifter(-1);
        pause(data.pause_between_sessions_s);
    end
    data.voice_onset_f = round(data.voice_onset_s * data.Fs / data.frameSize);
    data.voice_onset_ms = data.voice_onset_s * 1000;
    
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