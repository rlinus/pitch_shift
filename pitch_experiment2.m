rng('shuffle');

verbose = 1;

frameSize = 64;
Fs = 44100;

pitch_levels_cents = [-80 -40 40 80];
pitch_levels = 2.^(pitch_levels_cents/1200);
max_pitch_change=5;

num_sessions = 10;
session_duration_s = 1;
session_duration_f = round(session_duration_s * Fs / frameSize);
pitchchange_times_f = [0:session_duration_f:num_sessions*session_duration_f-1]';

pitch_level_sqs_cents = zeros(num_sessions, 1);
for i=1:num_sessions
    if rem(i,2) == 1
        pitch_level_sqs_cents(i) = 0;
    else
        lvl = ceil(length(pitch_levels_cents)*rand(1));
        pitch_level_sqs_cents(i) = pitch_levels_cents(lvl);
    end
end
        
pitch_level_sqs = 2.^(pitch_level_sqs_cents/1200);

rt_pitch_shifter(pitchchange_times_f,pitch_level_sqs);
pause(num_sessions*session_duration_s*1.01);
[y_r,y_ps] = rt_pitch_shifter(-1);

%%

y_r_session = cell(num_sessions,1);
mean_session = zeros(num_sessions,1);
for i=1:num_sessions
    y_r_session{i} = y_r((i-1)*session_duration_f*frameSize+1:i*session_duration_f*frameSize);
    [f0_time{i},f0_value{i}]=shrp(y_r_session{i},Fs);
    
    cand{i} = [false; abs(diff(f0_value{i}))<max_pitch_change] & f0_time{i} > 500;
      
    min_area = 3;
    c=0;
    for j=1:length(cand{i})
        if cand{i}(j) == true
            c = c+1;
            continue;
        else
            if c < min_area
                cand{i}(j-c:j-1) = false;
            end
            c=0;
        end
    end
    if c < min_area, cand{i}(end-c+1:end) = false; end;
    
    %apply Winsorized filter
    nout = 0.1;
    C = f0_time{i}(cand{i});
    [A,I] = sort(f0_value{i}(cand{i}));
    I_f = I(round(length(I)*nout/2):end-round(length(I)*nout/2));
    ca = false(length(cand{i}),1);
    for j=1:length(I_f)
        ca = ca | (f0_time{i} == C(I_f(j)));
    end
    cand{i} = cand{i} & ca;
    
    mean_session(i) = mean(f0_value{i}(cand{i}));
    std_session(i) = std(f0_value{i}(cand{i}));
    
    if verbose == 1
        figure(1);
        clf;
        
        t = linspace(0,length(y_r_session{i})*1000/Fs,length(y_r_session{i}));
        plot(t,y_r_session{i}*(max(f0_value{i})-mean_session(i))/max(y_r_session{i})+mean_session(i),'y');
        hold on
        
        plot(f0_time{i},f0_value{i}, 'b');
        hold on
        f0nan=nan(size(f0_value{i})); f0nan(cand{i})=f0_value{i}(cand{i});
        plot(f0_time{i},f0nan,'r');
        plot([min(f0_time{i}) max(f0_time{i})], [mean_session(i) mean_session(i)],'k');
        
        title(['Session ', num2str(i)]);
        xlabel('time [ms]');
        ylabel('frequency [Hz]');
        
        pause;
    end
end

pitch_changes = 1200 * log2(mean_session(2:end)./mean_session(1:end-1));

result = [pitch_level_sqs_cents(2:2:end), pitch_changes(1:2:end)]