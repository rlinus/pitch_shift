function data = vowel_exp_plot_mean(data)
    data.time_before_shift_ms = 300;
    data.time_after_shift_ms = 500;
    
    data.time_before_shift_ms = round(data.time_before_shift_ms/data.timestep)*data.timestep;
    data.time_after_shift_ms = round(data.time_after_shift_ms/data.timestep)*data.timestep;
    
    for i=1:data.num_sessions
        p =find(data.f0_time_s{i} >= data.voice_onset_ms(i)+data.shift_onset_ms(i),1);
        data.f0_value_s_cut{i} = data.f0_value_s{i}(p-(data.time_before_shift_ms/data.timestep):p+(data.time_after_shift_ms/data.timestep));
        data.voiced_regions_s_cut{i} = data.voiced_regions_s{i}(p-(data.time_before_shift_ms/data.timestep):p+(data.time_after_shift_ms/data.timestep));
        if ~all(data.voiced_regions_s_cut{i}), fprintf('Session %i is not completly voiced.\n', i); end;
        
        data.f0_value_ps_s_cut{i} = data.f0_value_ps_s{i}(p-(data.time_before_shift_ms/data.timestep):p+(data.time_after_shift_ms/data.timestep));

    end
    
    data.cols={'r','m','g','c','b','k','y'};
    
    data.f0_sums = zeros(length(data.pitch_levels_cents),length(data.f0_value_s_cut{1}));
    data.f0_qsums = data.f0_sums;
    data.f0_counts = data.f0_sums;
    
    data.f0_sums_ps = data.f0_sums;
    for i=1:data.num_sessions
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
    
    data.time = -data.time_before_shift_ms:data.timestep:data.time_after_shift_ms;
    h1 = figure;
    for i=1:length(data.pitch_levels_cents)
        subplot(3,1,1);
        plot(data.time,data.f0_means(i,:),'Color',data.cols{i},'LineWidth', 1);
        hold on;
        plot(data.time,data.f0_means_ps(i,:),'Color',data.cols{i},'LineWidth', 1,'LineStyle', '--');

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

end