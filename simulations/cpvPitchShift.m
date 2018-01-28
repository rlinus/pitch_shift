function [outdata,windowed_frame] = cpvPitchShift(pitchShift,input,init)
    %this is a translation of cpvPitchShift.cpp to MATLAB code for
    %simulation purposes

    if nargin < 3, init=false; end;
    
    fftFrameSize = 1024;
    stepSize = 64;
    sampleRate = 44100;
    
    osamp = fftFrameSize/stepSize;
    inFifoLatency = fftFrameSize-stepSize;
    
    unwrapdata = 2*pi*stepSize*(0:floor(fftFrameSize/2))/fftFrameSize;
    
    t=(1:fftFrameSize)-1;
    window = (-.5*cos(2*pi*t/fftFrameSize)+.5);
    
    input = input(:)';
    
    persistent gInFIFO;
    persistent gLastPhase;
    persistent gOutputAccum;
    persistent sphase;
    
    if isempty(gInFIFO) || init, gInFIFO=zeros(1,fftFrameSize);end;
    if isempty(gLastPhase) || init, gLastPhase=zeros(1,floor(fftFrameSize/2)+1);end;
    if isempty(gOutputAccum) || init, gOutputAccum=zeros(1,fftFrameSize*2);end;
    
%     SynthesisLen = round(pitchShift*stepSize);
%     Hopratio = SynthesisLen/stepSize;
%     ifftFrameSize=round(fftFrameSize/SynthesisLen*stepSize);
    ifftFrameSize = round(fftFrameSize/pitchShift);
    Hopratio = fftFrameSize/ifftFrameSize;
    
    
    gInFIFO(inFifoLatency+1:end) = input;
    
    gFFTworksp_real = gInFIFO.*window;
    
    gFFTworksp = fft(gFFTworksp_real);
    gFFTworksp = gFFTworksp(1:floor(fftFrameSize/2)+1);
    
    magn=abs(gFFTworksp);
    
    phase=angle(gFFTworksp);

    yunwrap = (phase - gLastPhase) - unwrapdata;
    gLastPhase = phase;
    yunwrap = yunwrap - round(yunwrap/(2*pi))*2*pi;
    
    yunwrap = (yunwrap + unwrapdata) * Hopratio;
    if init
        sphase = phase;
    else
        sphase = sphase + yunwrap;
    end

    gFFTworksp = magn.* complex(cos(sphase), sin(sphase));
    
    gFFTworksp_real = ifft(gFFTworksp,fftFrameSize,'symmetric');
    
    windowed_frame = gFFTworksp_real.*window/(osamp*3/8)*Hopratio;
    
    ywin=interpft(windowed_frame,ifftFrameSize);
    %%fprintf('l:%i\n',length(ywin));
    %%gain = 1/(fftFrameSize*sum(window.^2)/Hopratio/stepSize);
    
%     ywinc = zeros(1,fftFrameSize);
%     if(ifftFrameSize>=fftFrameSize)
%         b=floor((ifftFrameSize-fftFrameSize)/2);
%         ywinc = ywin(1+b:fftFrameSize+b);
%     else
%         b=floor((fftFrameSize-ifftFrameSize)/2);
%         ywinc(1+b:ifftFrameSize+b)=ywin;
%     end
%     gOutputAccum(1:fftFrameSize) = gOutputAccum(1:fftFrameSize)+ywinc;
    
    %%%%%gOutputAccum(ifftFrameSize-stepSize+1:end)=zeros(1,length(gOutputAccum(ifftFrameSize-stepSize+1:end)));
    
    gOutputAccum(1:ifftFrameSize) = gOutputAccum(1:ifftFrameSize)+ywin;
    
    outdata = gOutputAccum(1:stepSize);
    
    gOutputAccum = [gOutputAccum(stepSize+1:end), zeros(1,stepSize)];
    
    gInFIFO = [gInFIFO(stepSize+1:end), zeros(1,stepSize)];
end