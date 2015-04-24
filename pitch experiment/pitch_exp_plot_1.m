function data = pitch_exp_plot_1(data)
    data.ref_width = 20;

    data.num_shifts = floor(data.num_sessions_rec/2);

    data.refs = zeros(data.num_shifts,1);
    for i=1:data.num_shifts
        v = find(data.voiced_regions(data.f0_time<data.session_starts_ms(2*i)),data.ref_width,'last');
        data.refs(i) = mean(data.f0_value(v));
    end
    
    h = figure;
    plot(2:2:2*length(data.refs),data.refs);
    title('reference values [Hz]');
    
    f0nan=nan(size(data.f0_value)); f0nan(data.voiced_regions)=data.f0_value(data.voiced_regions);

    f = figure;
    p1 = subplot(7,1,[1 2 3]);
    p2 = subplot(7,1,[4 5 6]);

    sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',data.num_shifts,'Value',1,...
        'SliderStep', [1/(data.num_shifts-1) 1/(data.num_shifts-1)],...
        'Units','normalized',...
        'Position', [0.1 0.05 0.8 0.05],...
        'Callback', @plot_data);
    
    plot_data(sld);
    
    function plot_data(hObject, eventdata, handles)
        hObject.Value = round(hObject.Value);
        i = hObject.Value;
        
        t = linspace(-data.session_duration_ms,data.session_duration_ms,2*data.session_duration_f*data.frameSize+1);
        s = data.y_r(data.session_starts(2*i)-data.session_duration_f*data.frameSize:data.session_starts(2*i)+data.session_duration_f*data.frameSize);
        plot(p1,t,s*abs(diff(data.F0MinMax))/2/max(abs(s))+mean(data.F0MinMax),'y');
        
        hold(p1,'on');
        
        plot(p1,[-data.session_duration_ms, data.session_duration_ms],[data.refs(i),data.refs(i)],'k');
        
        if i~=data.num_shifts
            t_s = data.f0_time >= data.session_starts_ms(2*i-1) & data.f0_time < data.session_starts_ms(2*i+1);
        else
            t_s = data.f0_time >= data.session_starts_ms(2*i-1);
        end
        time = data.f0_time(t_s)-data.session_starts_ms(2*i);
        plot(p1,time,data.f0_value(t_s),'b');
        plot(p1,time,f0nan(t_s),'r');
        hold(p1,'off');
        
        ylim(p1,data.F0MinMax);
        xlim(p1,[-data.session_duration_ms, data.session_duration_ms]);
        ylabel(p1,'f0 [Hz]');
        
        plot(p2,time,1200*log2(f0nan(t_s)/data.refs(i)),'r');
        xlim(p2,[-data.session_duration_ms, data.session_duration_ms]);
        ylabel(p2,'pitch [cents]');
        xlabel(p2,'time [ms]');
        
        title(p1, sprintf('Session: %i, Applied Shift: %i cents',2*i,data.pitch_level_sqs_cents(2*i)));
    end
end
