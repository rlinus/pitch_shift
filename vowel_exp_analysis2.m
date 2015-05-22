function data = vowel_exp_analysis2(data)
    data.max_pitch_change = 10; %[Hz/ms]
    data.min_duration_of_voiced_regions = 20; %[ms]
    data.frame_length_f = 1024;
    data.timestep_f = 64;
    
    data.frame_length = 1000*data.frame_length_f/data.Fs; %[ms]
    data.timestep = 1000*data.timestep_f/data.Fs; %[ms]
    
    data.f0_low=0.75*data.piano_freq*2;
    data.f0_high=1.25*data.piano_freq*2; 


    
    for i=1:data.num_sessions
        [data.f0_time_s{i},data.f0_value_s{i}]=wavePitchOL(data.y_r{i},data.Fs,data.timestep_f,data.frame_length_f);
        [data.f0_time_ps_s{i},data.f0_value_ps_s{i}]=wavePitchOL(data.y_ps{i},data.Fs,data.timestep_f,data.frame_length_f);
        
%         pitch = data.detected_pitch{i}(1:length(data.f0_time_s{i}));
%         data.f0_value_s{i} = pitch;
%         data.f0_time_s{i} = data.f0_time_ps_s{i};
        
        for j=1:length(data.f0_value_s{i})
            if 2 * data.f0_value_s{i}(j) < data.f0_high && 2 * data.f0_value_s{i}(j) > data.f0_low
                data.f0_value_s{i}(j) = 2 * data.f0_value_s{i}(j);
            elseif 0.5 * data.f0_value_s{i}(j) < data.f0_high && 0.5 * data.f0_value_s{i}(j) > data.f0_low
                data.f0_value_s{i}(j) = 2 * data.f0_value_s{i}(j);
            end
        end
        
        
        data.voiced_regions_s{i} = [false; abs(diff(data.f0_value_s{i}))<data.max_pitch_change*data.timestep] & data.f0_value_s{i} < data.f0_high & data.f0_value_s{i} > data.f0_low;
        
        l=0;
        for j=1:length(data.voiced_regions_s{i})
            if data.voiced_regions_s{i}(j) == true
                l = l+1;
                continue;
            else
                if l < data.min_duration_of_voiced_regions/data.timestep;
                    data.voiced_regions_s{i}(j-l:j-1) = false;
                end
                l=0;
            end
        end
        if l < data.min_duration_of_voiced_regions/data.timestep, data.voiced_regions_s{i}(end-l+1:end) = false; end;
    end
end