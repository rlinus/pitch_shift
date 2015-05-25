function vowel_exp_plot_2c(data_c1,data_c2)
    
    if data_c1.piano_freq ~= data_c2.piano_freq 
        warning('reference frequencies must be the same'); return;
    end
    if data_c1.timestep ~= data_c2.timestep || data_c1.frame_length ~= data_c2.frame_length || data_c1.max_pitch_change ~= data_c2.max_pitch_change || data_c1.min_duration_of_voiced_regions ~= data_c2.min_duration_of_voiced_regions
        warning('analysis settings must be the same'); return;
    end

    if data_c1.mode ~= data_c2.mode
        warning('modes must be the same'); return;
    end
    
    if data_c1.time_before_shift_ms ~= data_c2.time_before_shift_ms || data_c1.time_after_shift_ms ~= data_c2.time_after_shift_ms
        warning('plot settings must be the same'); return;
    end
    
    cols={'r','m','g','c','b','k','y'};
    
    h1 = figure;
    h1a = subplot(2,1,1);
    h1b = subplot(2,1,2);
    
    for i=1:length(data_c1.pitch_levels_cents)
        subplot(2,1,1);
        hold on;
        
        plot(data_c1.time,data_c1.f0_means(i,:),'Color',cols{i},'LineWidth', 1);
        %plot(data_c1.time(l:end),data_c1.f0_means_ps(i,l:end),'Color',cols{i},'LineWidth', 1,'LineStyle', '--');

        subplot(2,1,2);
        plot(data_c1.time,data_c1.f0_std(i,:),'Color',cols{i},'LineWidth', 1);
        hold on;


        leg1{i} = sprintf('%i cents (%s)',data_c1.pitch_levels_cents(i), strrep(data_c1.condition_name, '_', ' '));
    end
    last_i = i;
    for i=1:length(data_c1.pitch_levels_cents)
        subplot(2,1,1);
        hold on;
        plot(data_c2.time,data_c2.f0_means(i,:),'Color',cols{i+last_i},'LineWidth', 1);
        %plot(data_c2.time(l:end),data_c2.f0_means_ps(i,l:end),'Color',cols{i+last_i},'LineWidth', 1,'LineStyle', '--');

        subplot(2,1,2);
        plot(data_c2.time,data_c2.f0_std(i,:),'Color',cols{i+last_i},'LineWidth', 1);
        hold on;

        leg1{i+last_i} = sprintf('%i cents (%s)',data_c2.pitch_levels_cents(i), strrep(data_c2.condition_name, '_', ' '));
    end
    xlim(h1a,[-data_c1.time_before_shift_ms,data_c1.time_after_shift_ms]);
    xlim(h1b,[-data_c1.time_before_shift_ms,data_c1.time_after_shift_ms]);
    ylabel(subplot(2,1,1),'f0 [Hz]'); ylabel(subplot(2,1,2),'std [Hz]');
    xlabel(subplot(2,1,1),'time [ms]'); xlabel(subplot(2,1,2),'time [ms]');
    legend(subplot(2,1,2),leg1);
    title(h1a,['subject: ', data_c1.subject]);
    
    plot(h1a,[0 0],h1a.YLim,'Color','k');
    plot(h1a,[data_c1.shift_duration_ms data_c1.shift_duration_ms],h1a.YLim,'Color','k');
    plot(h1b,[0 0],h1b.YLim,'Color','k');
    plot(h1b,[data_c1.shift_duration_ms data_c1.shift_duration_ms],h1b.YLim,'Color','k');
    
    h2 = figure;
    h2a = subplot(2,1,1);
    h2b = subplot(2,1,2);
    
    for i=1:length(data_c1.pitch_levels_cents)
        subplot(2,1,1);
        hold on;
        
        plot(data_c1.time,data_c1.f0_means_cents(i,:),'Color',cols{i},'LineWidth', 1);
        %plot(data_c1.time(l:end),data_c1.f0_means_ps_cents(i,l:end),'Color',cols{i},'LineWidth', 1,'LineStyle', '--');

        subplot(2,1,2);
        plot(data_c1.time,data_c1.f0_std(i,:),'Color',cols{i},'LineWidth', 1);
        hold on;


        leg1{i} = sprintf('%i cents (%s)',data_c1.pitch_levels_cents(i), strrep(data_c1.condition_name, '_', ' '));
    end
    last_i = i;
    for i=1:length(data_c1.pitch_levels_cents)
        subplot(2,1,1);
        hold on;
        plot(data_c2.time,data_c2.f0_means_cents(i,:),'Color',cols{i+last_i},'LineWidth', 1);
        %plot(data_c2.time(l:end),data_c2.f0_means_ps_cents(i,l:end),'Color',cols{i+last_i},'LineWidth', 1,'LineStyle', '--');

        subplot(2,1,2);
        plot(data_c2.time,data_c2.f0_std(i,:),'Color',cols{i+last_i},'LineWidth', 1);
        hold on;

        leg1{i+last_i} = sprintf('%i cents (%s)',data_c2.pitch_levels_cents(i), strrep(data_c2.condition_name, '_', ' '));
    end
    xlim(h2a,[-data_c1.time_before_shift_ms,data_c1.time_after_shift_ms]);
    xlim(h2b,[-data_c1.time_before_shift_ms,data_c1.time_after_shift_ms]);
    ylabel(subplot(2,1,1),'f0 [cents]'); ylabel(subplot(2,1,2),'std [Hz]');
    xlabel(subplot(2,1,1),'time [ms]'); xlabel(subplot(2,1,2),'time [ms]');
    legend(subplot(2,1,2),leg1);
    title(h2a,['subject: ', data_c1.subject]);
    
    plot(h2a,[0 0],h2a.YLim,'Color','k');
    plot(h2a,[data_c1.shift_duration_ms data_c1.shift_duration_ms],h2a.YLim,'Color','k');
    plot(h2b,[0 0],h2b.YLim,'Color','k');
    plot(h2b,[data_c1.shift_duration_ms data_c1.shift_duration_ms],h2b.YLim,'Color','k');


end

