function data = vowel_exp_plot(data)
    data.correct_bias = true;
    data.ref = 'constant';
    %data.ref = 'time variant';
    data.invalid_sessions = [];
    data.ref_freq = data.piano_freq;
    data.latency_ms=20;
    
    data.time_before_shift_ms = 300;
    data.time_after_shift_ms = 1300;
    
    if data.mode==1, data.time_before_shift_ms=0; end;

    data.time_before_shift_ms = round(data.time_before_shift_ms/data.timestep)*data.timestep;
    data.time_after_shift_ms = round(data.time_after_shift_ms/data.timestep)*data.timestep;
    before = round(data.time_before_shift_ms/data.timestep);
    after =  round(data.time_after_shift_ms/data.timestep);
    data.time = -data.time_before_shift_ms:data.timestep:data.time_after_shift_ms;
    
    data.valid_sessions = 1:data.num_sessions;
    data.valid_sessions(data.invalid_sessions) = [];
    
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
    hold(p1,'on'); hold(p2,'on');
    
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
        t = (1:length(data.static_pitch_factor_sqs{i})) * data.frameSize * 1000 / data.Fs;
        
        
        cla(p1);
        ylim(p1,[data.f0_low data.f0_high]);
        xlim(p1,[0, data.f0_time_s{i}(end)]);
        
        plot(p1,p1.XLim,[data.ref_freq,data.ref_freq],'k','LineStyle', '--');
        
        plot(p1,[data.voice_onset_ms(i), data.voice_onset_ms(i)],p1.YLim,'k');
        
        if(data.mode==2)
            plot(p1,[data.voice_onset_ms(i)+data.shift_onset_ms(i), data.voice_onset_ms(i)+data.shift_onset_ms(i)],p1.YLim,'k','LineStyle', '--');
            plot(p1,[data.voice_onset_ms(i)+data.shift_onset_ms(i)-data.time_before_shift_ms, data.voice_onset_ms(i)+data.shift_onset_ms(i)-data.time_before_shift_ms],p1.YLim,'y');
            plot(p1,[data.voice_onset_ms(i)+data.shift_onset_ms(i)+data.time_after_shift_ms, data.voice_onset_ms(i)+data.shift_onset_ms(i)+data.time_after_shift_ms],p1.YLim,'y');
        end

        plot(p1,t,data.var_pitch_factor_sqs{i}.*data.static_pitch_factor_sqs{i}.*data.control_pitch_factor_sqs{i}*data.ref_freq,'Color',[0.6 0.6 0.6]);
        plot(p1,t,data.static_pitch_factor_sqs{i}.*data.control_pitch_factor_sqs{i}*data.ref_freq,'k');
        
        plot(p1,t,data.detected_pitch{i},'c');
        
        f0nan=nan(size(data.f0_value_s{i})); f0nan(data.voiced_regions_s{i})=data.f0_value_s{i}(data.voiced_regions_s{i});
        plot(p1,data.f0_time_s{i},data.f0_value_s{i},'b');
        plot(p1,data.f0_time_s{i},f0nan,'r','LineWidth', 1.5);
        plot(p1,data.f0_time_s{i},data.f0_value_ps_s{i},'g');

        
        
        
        ylabel(p1,'f0 [Hz]');
        
        cla(p2);
        ylim(p2,[-1*max(abs(data.pitch_levels_cents))-150, max(abs(data.pitch_levels_cents))+150]);
        xlim(p2,[0, data.f0_time_s{i}(end)]);
        
        plot(p2,p2.XLim,[0,0],'k','LineStyle', '--');
        
        plot(p2,[data.voice_onset_ms(i), data.voice_onset_ms(i)],p2.YLim,'k');
        if(data.mode==2)
            plot(p2,[data.voice_onset_ms(i)+data.shift_onset_ms(i), data.voice_onset_ms(i)+data.shift_onset_ms(i)],p2.YLim,'k','LineStyle', '--');
        end
        
        plot(p2,t,1200*log2(data.var_pitch_factor_sqs{i}.*data.static_pitch_factor_sqs{i}.*data.control_pitch_factor_sqs{i}),'Color',[0.6 0.6 0.6]);
        plot(p2,t,1200*log2(data.static_pitch_factor_sqs{i}.*data.control_pitch_factor_sqs{i}),'k');
        
        plot(p2,t,1200*log2(data.detected_pitch{i}/data.ref_freq),'c');
        
        plot(p2,data.f0_time_s{i},1200*log2(data.f0_value_s{i}/data.ref_freq),'b');
        plot(p2,data.f0_time_s{i},1200*log2(f0nan/data.ref_freq),'r','LineWidth', 1);
        plot(p2,data.f0_time_s{i},1200*log2(data.f0_value_ps_s{i}/data.ref_freq),'g');

        
        
        ylabel(p2,'pitch [cents]');
        xlabel(p2,'time [ms]');
        
        title(p1, sprintf('Session: %i, Applied Shift: %i cents',i,data.pitch_level_sqs_cents(i)));
    end

    %% plot mean over sessions
    
    cols={'r','m','g','c','b','k','y'};
    l = round(data.latency_ms/data.timestep)+1;
    
    data.f0_sums = zeros(length(data.pitch_levels_cents),length(data.f0_value_s_cut{1}));
    data.f0_qsums = data.f0_sums;
    data.f0_counts = data.f0_sums;
    
    data.f0_sums_ps = data.f0_sums;
    
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
    
    if data.correct_bias == true
        for i=1:length(data.pitch_levels_cents)
            data.f0_means(i,:)= data.f0_means(i,:)/mean(data.f0_means(i,1:before))*data.ref_freq;
            data.f0_means_ps(i,:)= data.f0_means_ps(i,:)/mean(data.f0_means_ps(i,1:before))*data.ref_freq;
        end
    end
    
    h1 = figure;
    for i=1:length(data.pitch_levels_cents)
        subplot(3,1,1);
        plot(data.time,data.f0_means(i,:),'Color',cols{i},'LineWidth', 1);
        hold on;
        plot(data.time(l:end),data.f0_means_ps(i,l:end),'Color',cols{i},'LineWidth', 1,'LineStyle', '--');

        subplot(3,1,2);
        plot(data.time,data.f0_std(i,:),'Color',cols{i},'LineWidth', 1);
        hold on;

        subplot(3,1,3);
        plot(data.time,data.f0_counts(i,:),'Color',cols{i},'LineWidth', 1);
        hold on;

        leg1{i} = sprintf('%i cents',data.pitch_levels_cents(i));
    end
    xlim(subplot(3,1,1),[-data.time_before_shift_ms,data.time_after_shift_ms]);
    xlim(subplot(3,1,2),[-data.time_before_shift_ms,data.time_after_shift_ms]);
    xlim(subplot(3,1,3),[-data.time_before_shift_ms,data.time_after_shift_ms]);
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
        plot(data.time,data.f0_means_cents(i,:),'Color',cols{i},'LineWidth', 1);
        hold on;
        plot(data.time(l:end),data.f0_means_ps_cents(i,l:end),'Color',cols{i},'LineWidth', 1,'LineStyle', '--');

        subplot(3,1,2);
        plot(data.time,data.f0_std(i,:),'Color',cols{i},'LineWidth', 1);
        hold on;

        subplot(3,1,3);
        plot(data.time,data.f0_counts(i,:),'Color',cols{i},'LineWidth', 1);
        hold on;

        leg2{end+1} = sprintf('%i cents',data.pitch_levels_cents(i));
    end
    xlim(subplot(3,1,1),[-data.time_before_shift_ms,data.time_after_shift_ms]);
    xlim(subplot(3,1,2),[-data.time_before_shift_ms,data.time_after_shift_ms]);
    xlim(subplot(3,1,3),[-data.time_before_shift_ms,data.time_after_shift_ms]);
    ylabel(subplot(3,1,1),'f0 [cents]'); ylabel(subplot(3,1,2),'std [Hz]'); ylabel(subplot(3,1,3),'number of trials');
    xlabel(subplot(3,1,1),'time [ms]'); xlabel(subplot(3,1,2),'time [ms]'); xlabel(subplot(3,1,3),'time [ms]');
    %ylim(subplot(3,1,1),[2*min(data.pitch_levels_cents) 2*max(data.pitch_levels_cents)]);
    legend(subplot(3,1,3),leg2);
end

