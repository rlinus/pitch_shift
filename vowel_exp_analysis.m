function data = vowel_exp_analysis(data)
    data.sex = 'female';
    data.max_pitch_change = 0.5; %[Hz/ms]
    data.frame_length = 2.5*40; %[ms]
    data.timestep = 10/2; %[ms]
    data.min_duration_of_voiced_regions = 150; %[ms]

    if strcmpi(data.sex,'male')
        data.F0MinMax = [50 250];
        data.f0_low=100;
        data.f0_high=150; 
    else
        data.F0MinMax = [300 500];
        data.f0_low=350;
        data.f0_high=450;
    end
    
    for i=1:data.num_sessions
        [data.f0_time_s{i},data.f0_value_s{i}]=shrp(data.y_r{i},data.Fs,data.F0MinMax,data.frame_length,data.timestep);
        [data.f0_time_ps_s{i},data.f0_value_ps_s{i}]=shrp(data.y_ps{i},data.Fs,data.F0MinMax,data.frame_length,data.timestep);
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