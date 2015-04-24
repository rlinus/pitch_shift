    function data = pitch_exp_plot_2(data)
    data.analysis_onset = 'first_voiced_region'; %or 'first_sample'
    data.cols={'r','m','g','c','b','k','y'};

    data.pitch_levels_cents_wz = [data.pitch_levels_cents(:); 0];

    data.f0_time_s = cell(data.num_sessions_rec,1);
    data.f0_value_s = cell(data.num_sessions_rec,1);
    data.voiced_regions_s = cell(data.num_sessions_rec,1);

    for i=1:data.num_sessions_rec-1
        data.y_r_s{i} = data.y_r(data.session_starts(i):data.session_starts(i+1)-1);
        data.y_ps_s{i} = data.y_ps(data.session_starts(i):data.session_starts(i+1)-1);

        t_s = data.f0_time >= data.session_starts_ms(i) & data.f0_time < data.session_starts_ms(i+1);
        data.f0_time_s{i} = data.f0_time(t_s);
        data.f0_value_s{i} = data.f0_value(t_s);
        data.voiced_regions_s{i} = data.voiced_regions(t_s);
    end
    data.y_r_s{data.num_sessions_rec} = data.y_r(data.session_starts(end):end);
    data.y_ps_s{data.num_sessions_rec} = data.y_ps_s(data.session_starts(end):end);
    t_s = data.f0_time >= data.session_starts_ms(end);
    data.f0_time_s{data.num_sessions_rec} = data.f0_time(t_s);
    data.f0_value_s{data.num_sessions_rec} = data.f0_value(t_s);
    data.voiced_regions_s{data.num_sessions_rec} = data.voiced_regions(t_s);

    data.t_on = zeros(data.num_sessions_rec,1);
    for i=1:data.num_sessions_rec
        if strcmpi(data.analysis_onset,'first_voiced_region')
            try
                data.t_on(i) = find(data.voiced_regions_s{i},1);
            catch
                data.t_on(i) = 1;
            end
        else
            data.t_on(i) = 1;
        end

        %data.len(i) = length(data.voiced_regions_s{i}(data.t_on(i):end));
    end

    %data.f0_sums = zeros(length(data.pitch_levels_cents_wz),min(data.len));
    data.f0_sums = zeros(length(data.pitch_levels_cents_wz),floor(data.session_duration_ms/data.timestep));
    data.f0_qsums = data.f0_sums;
    data.f0_counts = data.f0_sums;
    for i=1:data.num_sessions_rec
        pos = find(data.pitch_levels_cents_wz == data.pitch_level_sqs_cents(i),1);

        for j=1:size(data.f0_sums,2);
            if length(data.voiced_regions_s{i}) >= data.t_on(i)-1+j && data.voiced_regions_s{i}(data.t_on(i)-1+j)
                data.f0_sums(pos,j) = data.f0_sums(pos,j) + data.f0_value_s{i}(data.t_on(i)-1+j);
                data.f0_qsums(pos,j) = data.f0_qsums(pos,j) + data.f0_value_s{i}(data.t_on(i)-1+j).^2;
                data.f0_counts(pos,j) = data.f0_counts(pos,j) +1;
            end
        end
    end

    data.f0_means = data.f0_sums ./ data.f0_counts;
    data.f0_std = sqrt(data.f0_qsums ./ data.f0_counts - data.f0_means.^2);


    refs = repmat(data.f0_means(end,:),length(data.pitch_levels_cents),1);
    data.pitches = 1200 * log2(data.f0_means(1:end-1,:)./refs);

    data.pitches_percent_compensation = data.pitches ./ repmat(-1*data.pitch_levels_cents(:),1,size(data.pitches,2))*100;


    data.x = 0:data.timestep:(size(data.f0_sums,2)-1)*data.timestep;
    h1 = figure;
    for i=1:length(data.pitch_levels_cents_wz)
        subplot(3,1,1);
        plot(data.x,data.f0_means(i,:),'Color',data.cols{i},'LineWidth', 1);
        hold on;

        subplot(3,1,2);
        plot(data.x,data.f0_std(i,:),'Color',data.cols{i},'LineWidth', 1);
        hold on;

        subplot(3,1,3);
        plot(data.x,data.f0_counts(i,:),'Color',data.cols{i},'LineWidth', 1);
        hold on;

        leg1{i} = sprintf('%i cents',data.pitch_levels_cents_wz(i));
    end
    ylabel(subplot(3,1,1),'f0 [Hz]'); ylabel(subplot(3,1,2),'std [Hz]'); ylabel(subplot(3,1,3),'number of trials');
    xlabel(subplot(3,1,1),'time [ms]'); xlabel(subplot(3,1,2),'time [ms]'); xlabel(subplot(3,1,3),'time [ms]');
    legend(subplot(3,1,3),leg1);

    h2 = figure;
    for i=1:length(data.pitch_levels_cents)
        subplot(2,1,1);
        hold on;
        plot(data.x,data.pitches(i,:)+data.pitch_levels_cents(i),'Color',data.cols{i},'LineWidth', 1);
        %plot([x(1), x(end)],[-data.pitch_levels_cents(i), -data.pitch_levels_cents(i)],'Color',cols{i},'LineStyle','--');

        subplot(2,1,2);
        hold on;
        plot(data.x,data.pitches_percent_compensation(i,:),'Color',data.cols{i},'LineWidth', 1);

        leg2{i} = sprintf('%i cents',data.pitch_levels_cents(i));
    end
    ylabel(subplot(2,1,1),'pitch shift [cents]'); ylabel(subplot(2,1,2),'relative pitch shift [%]');
    xlabel(subplot(2,1,1),'time [ms]'); xlabel(subplot(2,1,2),'time [ms]');
    legend(subplot(2,1,2),leg2);
end