function [perceived_shift, invalid_trial] = get_perceived_shift(signal, init_shift)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    step_size = 10;

    sel_shift = init_shift;
    
    pitch_loop(1,signal,2^(sel_shift/1200),2,0);

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
    pitch_loop(-1);
    
    perceived_shift = sel_shift;
    
    

    function wheel_clb(src,callbackdata)
        sel_shift = sel_shift - callbackdata.VerticalScrollCount * step_size;
        pitch_loop(0,2^(sel_shift/1200));
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

