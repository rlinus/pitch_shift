function [perceived_shift, invalid_trial] = get_perceived_shift(signal, init_shift, params)
%get_perceived_shift plays signal in endless loop and opens gui in which
%the subject can dynamically choose a pitch shift with the mouse wheel.

    step_size = 10; %increment size of the pitch shift

    sel_shift = init_shift;
    params.pitch_factor = init_shift;
    
    PitchLoop(1,params, signal);

    invalid_trial=0;

    ui = figure;

    ui.WindowScrollWheelFcn = @wheel_clb;

    ui.MenuBar = 'none';

    msg = {'Use the mouse wheel to set the pitch and press OK afterwards.',...
        ' ',...
        'Up: higher pitch',...
        'Down: lower pitch'};

    txt = uicontrol(ui, 'Style','text',...
            'Units', 'normalized',...
            'Position',[0 0.75 1 0.25],...
            'String',msg,...
            'FontSize',14);

    btn = uicontrol('Style', 'pushbutton', 'String', 'OK',...
        'Units', 'normalized',...
        'Position', [0.55 0.1 0.2 0.1],...
        'Callback', @ok_clb);
    
    btn = uicontrol('Style', 'pushbutton', 'String', 'Invalid Trial',...
        'Units', 'normalized',...
        'Position', [0.25 0.1 0.2 0.1],...
        'Callback', @invalid_clb);


    uiwait(ui);
    PitchLoop(-1);
    
    perceived_shift = sel_shift;
    
    

    function wheel_clb(src,callbackdata)
        sel_shift = sel_shift - callbackdata.VerticalScrollCount * step_size;
        PitchLoop(0,2^(sel_shift/1200));
    end

    function ok_clb(src,callbackdata)
        uiresume(ui);
        close(ui);
    end

    function invalid_clb(src,callbackdata)
        invalid_trial = 1;
        uiresume(ui);
        close(ui);
    end

end

