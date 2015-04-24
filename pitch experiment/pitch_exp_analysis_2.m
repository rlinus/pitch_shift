function data = pitch_exp_analysis_2(data)
    threshold = 0.5;

    data.F0MinMax = [50 250];
    
    sig = Signal(data.y_r', data.Fs);
    
    data.gs = sig.getOnsets(50,250);
    %r = sig.getOnsets(data.F0MinMax(1),data.F0MinMax(2));
    data.voiced_regions = r>threshold;
    
    [data.f0_value]= sig.mainPitch(data.F0MinMax(1), data.F0MinMax(2));
    data.f0_time = sig.framesPositions *1000/data.Fs;
end