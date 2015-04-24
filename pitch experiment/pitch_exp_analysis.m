function data = pitch_exp_analysis(data)
    data.sex = 'male';
    data.max_pitch_change = 4;
    data.frame_length = 40; %[ms]
    data.timestep = 10; %[ms]
    data.min_duration_of_voiced_regions = 150; %[ms]

    if strcmpi(data.sex,'male')
        data.F0MinMax = [50 250];
        data.f0_low=100;
        data.f0_high=150; 
    else
        data.F0MinMax = [120 400];
        data.f0_low=120;
        data.f0_high=350;
    end

    data.y_r_s = cell(data.num_sessions_rec,1);
    data.y_ps_s = data.y_r_s;


    [data.f0_time,data.f0_value]=shrp(data.y_r,data.Fs,data.F0MinMax,data.frame_length,data.timestep);

    data.voiced_regions = [false; abs(diff(data.f0_value))<data.max_pitch_change] & data.f0_value < data.f0_high & data.f0_value > data.f0_low;

    l=0;
    for j=1:length(data.voiced_regions)
        if data.voiced_regions(j) == true
            l = l+1;
            continue;
        else
            if l < data.min_duration_of_voiced_regions/data.timestep;
                data.voiced_regions(j-l:j-1) = false;
            end
            l=0;
        end
    end
    if l < data.min_duration_of_voiced_regions/data.timestep, data.voiced_regions(end-l+1:end) = false; end;
end

