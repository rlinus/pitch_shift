function data = vowel_exp_plot(data)

    data.ref = 'constant';
    %data.ref = 'time variant';
    data.invalid_sessions = [];
    data.ref_freq = 125;%data.piano_freq;
    data.latency_ms=20;
    if data.mode==1
        data.time_after_voice_onset_ms = 2800;
        
        data.time_after_voice_onset_ms = round(data.time_after_voice_onset_ms/data.timestep)*data.timestep;
        before = 0;
        after = data.time_after_voice_onset_ms/data.timestep;
        data.time = 0:data.timestep:data.time_after_voice_onset_ms;
    else
        data.time_before_shift_ms = 200;
        data.time_after_shift_ms = 1300;
        
        data.time_before_shift_ms = round(data.time_before_shift_ms/data.timestep)*data.timestep;
        data.time_after_shift_ms = round(data.time_after_shift_ms/data.timestep)*data.timestep;
        before = data.time_before_shift_ms/data.timestep;
        after =  data.time_after_shift_ms/data.timestep;
        data.time = -data.time_before_shift_ms:data.timestep:data.time_after_shift_ms;
    end
    
    for i=1:data.num_sessions
        if(data.mode==1) %shift before voice onset
            p =find(data.voiced_regions_s{i} & (data.f0_time_s{i}>=data.voice_onset_ms(i)),1);
            if isempty(p)
                data.invalid_sessions(end+1) = i;
                fprintf('session %i is invalid: no voiced regions.\n',i);
                continue;
            end;
        else  %shift after voice onset
            p =find(data.f0_time_s{i} >= data.voice_onset_ms(i)+data.shift_onset_ms(i),1);   
        end
        
        if(p-before < 1 || p+after > length(data.f0_value_s{i}))
            data.invalid_sessions(end+1) = i;
            fprintf('session %i is invalid: time window is out of recorded range.\n',i);
            continue;
        end
        data.f0_value_s_cut{i} = data.f0_value_s{i}(p-before:p+after);
        data.voiced_regions_s_cut{i} = data.voiced_regions_s{i}(p-before:p+after);
        if ~all(data.voiced_regions_s_cut{i}), fprintf('Session %i is not completly voiced.\n', i); end;

        data.f0_value_ps_s_cut{i} = data.f0_value_ps_s{i}(p-before:p+after);
    end
        

    %% plot seperate sessions
    f = figure;
    p1 = subplot(7,1,[1 2 3]);
    p2 = subplot(7,1,[4 5 6]);
    
    if data.num_sessions>1
        sld = uicontrol('Style', 'slider',...
            'Min',1,'Max',data.num_sessions,'Value',1,...
            'SliderStep', [1/(data.num_sessions-1) 1/(data.num_sessions-1)],...
            'Units','normalized',...
            'Position', [0.1 0.05 0.8 0.05],...
            'Callback', @plot_data);
        
        plot_data(sld);
    else
        plot_data
    end
    
    
    
    function plot_data(hObject, eventdata, handles)
        if nargin > 0
            hObject.Value = round(hObject.Value);
            i = hObject.Value;
        else
            i=1;
        end
        
%         t = linspace(-data.session_duration_ms,data.session_duration_ms,2*data.session_duration_f*data.frameSize+1);
%         s = data.y_r(data.session_starts(2*i)-data.session_duration_f*data.frameSize:data.session_starts(2*i)+data.session_duration_f*data.frameSize);
%         plot(p1,t,s*abs(diff(data.F0MinMax))/2/max(abs(s))+mean(data.F0MinMax),'y');
%         
%         hold(p1,'on');
        
        plot(p1,[data.f0_time_s{i}(1), data.f0_time_s{i}(end)],[data.ref_freq,data.ref_freq],'k');
        hold(p1,'on');
        plot(p1,[data.voice_onset_ms(i), data.voice_onset_ms(i)],[0.9*data.ref_freq 1.1*data.ref_freq],'k');
        if(data.mode==2)
            plot(p1,[data.voice_onset_ms(i)+data.shift_onset_ms(i), data.voice_onset_ms(i)+data.shift_onset_ms(i)],[0.9*data.ref_freq 1.1*data.ref_freq],'k','LineStyle', '--');
        end
%         if i~=data.num_shifts
%             t_s = data.f0_time >= data.session_starts_ms(2*i-1) & data.f0_time < data.session_starts_ms(2*i+1);
%         else
%             t_s = data.f0_time >= data.session_starts_ms(2*i-1);
%         end
%         time = data.f0_time(t_s)-data.session_starts_ms(2*i);
        f0nan=nan(size(data.f0_value_s{i})); f0nan(data.voiced_regions_s{i})=data.f0_value_s{i}(data.voiced_regions_s{i});
        plot(p1,data.f0_time_s{i},data.f0_value_s{i},'b');
        plot(p1,data.f0_time_s{i},f0nan,'r');
        plot(p1,data.f0_time_s{i},data.f0_value_ps_s{i},'g');
        hold(p1,'off');
        
        ylim(p1,[0.9*data.ref_freq 1.1*data.ref_freq]);
        xlim(p1,[data.f0_time_s{i}(1), data.f0_time_s{i}(end)]);
        ylabel(p1,'f0 [Hz]');
        
        plot(p2,data.f0_time_s{i},1200*log2(data.f0_value_s{i}/data.ref_freq),'b');
        hold(p2,'on');
        plot(p2,data.f0_time_s{i},1200*log2(data.f0_value_ps_s{i}/data.ref_freq),'g');
        hold(p2,'off');
        ylim(p2,[-2*max([abs(data.pitch_levels_cents),50]) 2*max([abs(data.pitch_levels_cents),50])]);
        xlim(p2,[data.f0_time_s{i}(1), data.f0_time_s{i}(end)]);
        ylabel(p2,'pitch [cents]');
        xlabel(p2,'time [ms]');
        
        title(p1, sprintf('Session: %i, Applied Shift: %i cents',i,data.pitch_level_sqs_cents(i)));
    end

    %% plot mean over sessions
    
    data.cols={'r','m','g','c','b','k','y'};
    l = round(data.latency_ms/data.timestep)+1;
    
    data.f0_sums = zeros(length(data.pitch_levels_cents),length(data.f0_value_s_cut{1}));
    data.f0_qsums = data.f0_sums;
    data.f0_counts = data.f0_sums;
    
    data.f0_sums_ps = data.f0_sums;
    
    data.valid_sessions = 1:data.num_sessions;
    data.valid_sessions(data.invalid_sessions) = [];
    for i=data.valid_sessions
        pos = find(data.pitch_levels_cents == data.pitch_level_sqs_cents(i),1);

        for j=1:size(data.f0_sums,2);
            if data.voiced_regions_s_cut{i}(j)
                data.f0_sums(pos,j) = data.f0_sums(pos,j) + data.f0_value_s_cut{i}(j);
                data.f0_qsums(pos,j) = data.f0_qsums(pos,j) + data.f0_value_s_cut{i}(j).^2;
                data.f0_counts(pos,j) = data.f0_counts(pos,j) + 1;
                
                
                data.f0_sums_ps(pos,j) = data.f0_sums_ps(pos,j) + data.f0_value_ps_s_cut{i}(j);
            end
        end
    end

    data.f0_means = data.f0_sums ./ data.f0_counts;
    data.f0_std = sqrt(data.f0_qsums ./ data.f0_counts - data.f0_means.^2);
    
    data.f0_means_ps = data.f0_sums_ps ./ data.f0_counts;
    
    h1 = figure;
    for i=1:length(data.pitch_levels_cents)
        subplot(3,1,1);
        plot(data.time,data.f0_means(i,:),'Color',data.cols{i},'LineWidth', 1);
        hold on;
        plot(data.time(l:end),data.f0_means_ps(i,l:end),'Color',data.cols{i},'LineWidth', 1,'LineStyle', '--');

        subplot(3,1,2);
        plot(data.time,data.f0_std(i,:),'Color',data.cols{i},'LineWidth', 1);
        hold on;

        subplot(3,1,3);
        plot(data.time,data.f0_counts(i,:),'Color',data.cols{i},'LineWidth', 1);
        hold on;

        leg1{i} = sprintf('%i cents',data.pitch_levels_cents(i));
    end
    ylabel(subplot(3,1,1),'f0 [Hz]'); ylabel(subplot(3,1,2),'std [Hz]'); ylabel(subplot(3,1,3),'number of trials');
    xlabel(subplot(3,1,1),'time [ms]'); xlabel(subplot(3,1,2),'time [ms]'); xlabel(subplot(3,1,3),'time [ms]');
    legend(subplot(3,1,3),leg1);
    
    
    if strcmpi(data.ref, 'time variant')
        ref_values = data.f0_means(find(data.pitch_levels_cents == 0,1),:);
        range = find(data.pitch_levels_cents ~= 0);
    else
        ref_values = data.ref_freq*ones(1,size(data.f0_sums,2));
        range = 1:length(data.pitch_levels_cents);
    end
    
    data.f0_means_cents = 1200 * log2(data.f0_means./repmat(ref_values,size(data.f0_means,1),1));
    data.f0_means_ps_cents = 1200 * log2(data.f0_means_ps./repmat(ref_values,size(data.f0_means_ps,1),1));
    
    leg2=cell(0);
    h2 = figure;
    for i=range
        subplot(3,1,1);
        plot(data.time,data.f0_means_cents(i,:),'Color',data.cols{i},'LineWidth', 1);
        hold on;
        plot(data.time(l:end),data.f0_means_ps_cents(i,l:end),'Color',data.cols{i},'LineWidth', 1,'LineStyle', '--');

        subplot(3,1,2);
        plot(data.time,data.f0_std(i,:),'Color',data.cols{i},'LineWidth', 1);
        hold on;

        subplot(3,1,3);
        plot(data.time,data.f0_counts(i,:),'Color',data.cols{i},'LineWidth', 1);
        hold on;

        leg2{end+1} = sprintf('%i cents',data.pitch_levels_cents(i));
    end
    ylabel(subplot(3,1,1),'f0 [cents]'); ylabel(subplot(3,1,2),'std [Hz]'); ylabel(subplot(3,1,3),'number of trials');
    xlabel(subplot(3,1,1),'time [ms]'); xlabel(subplot(3,1,2),'time [ms]'); xlabel(subplot(3,1,3),'time [ms]');
    %ylim(subplot(3,1,1),[2*min(data.pitch_levels_cents) 2*max(data.pitch_levels_cents)]);
    legend(subplot(3,1,3),leg2);
end

