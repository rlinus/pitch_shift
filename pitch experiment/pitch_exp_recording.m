function data = pitch_exp_recording
    rng('shuffle');

    data.frameSize = 64;
    data.Fs = 44100;

    data.pitch_levels_cents = [-60 -30 30 60];
    data.pitch_levels = 2.^(data.pitch_levels_cents/1200);

    data.num_sessions = 10;
    data.session_duration_ms = 1000;
    data.low_amp_threshold = 0.01;
    data.min_low_amp_duration_ms = 5;


    data.min_low_amp_duration_f = round(data.min_low_amp_duration_ms*data.Fs/data.frameSize / 1000);
    data.session_duration_f = round(data.session_duration_ms * data.Fs / data.frameSize / 1000);


    data.pitch_level_sqs_cents = zeros(data.num_sessions, 1);
    for i=1:data.num_sessions
        if rem(i,2) == 1
            data.pitch_level_sqs_cents(i) = 0;
        else
            lvl = ceil(length(data.pitch_levels_cents)*rand(1));
            data.pitch_level_sqs_cents(i) = data.pitch_levels_cents(lvl);
        end
    end

    data.pitch_level_sqs = 2.^(data.pitch_level_sqs_cents/1200);

    rt_pitch_shifter(data.session_duration_f,data.pitch_level_sqs,data.min_low_amp_duration_f ,data.low_amp_threshold,0,0.002);
    while(rt_pitch_shifter(0)==0)
        pause(0.5);
    end
    last_sample = rt_pitch_shifter(0);
    [data.y_r,data.y_ps,data.session_starts] = rt_pitch_shifter(-1);

    data.y_r = data.y_r(1:last_sample);
    data.y_ps = data.y_ps(1:last_sample);

    data.session_starts_ms = (data.session_starts-1)*1000/data.Fs;
    data.num_sessions_rec = length(data.session_starts);
end