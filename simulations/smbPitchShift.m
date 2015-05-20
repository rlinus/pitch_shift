function [outdata,windowed_frame] = smbPitchShift(pitchShift,input,init)
    if nargin < 3, init=false; end;
    
    fftFrameSize = 1024;
    stepSize =64;
    sampleRate = 44100;
    
    osamp = fftFrameSize/stepSize;
    inFifoLatency = fftFrameSize-stepSize;
    
    freqPerBin = sampleRate/fftFrameSize;
    expct = 2*pi*stepSize/fftFrameSize;
    
    persistent gInFIFO;
    persistent gLastPhase;
    persistent gSumPhase;
    persistent gOutputAccum;
    persistent lastShift;
    
    t=(1:fftFrameSize)-1;
    window = (-.5*cos(2*pi*t/fftFrameSize)+.5);

    if isempty(gInFIFO) || init, gInFIFO=zeros(1,fftFrameSize);end;
    if isempty(gLastPhase) || init, gLastPhase=zeros(1,floor(fftFrameSize/2)+1);end;
    if isempty(gSumPhase) || init, gSumPhase=zeros(1,floor(fftFrameSize/2)+1);end;
    if isempty(gOutputAccum) || init, gOutputAccum=zeros(1,fftFrameSize);end;
    if isempty(lastShift) || init, lastShift=0;end;
    
%     if(lastShift~=pitchShift)
%         %gSumPhase=gLastPhase;
%         %gSumPhase=zeros(1,floor(fftFrameSize/2)+1);
%         gLastPhase=zeros(1,floor(fftFrameSize/2)+1);
%                 
%     end;
    
    input = input(:)';
    
    gInFIFO(inFifoLatency+1:end) = input;
    
    gFFTworksp_real = gInFIFO.*window;
    
    gFFTworksp = fft(gFFTworksp_real);
    gFFTworksp = gFFTworksp(1:floor(fftFrameSize/2)+1);
    
    
    magn=abs(gFFTworksp);
    
    phase=angle(gFFTworksp);
    
    tmp = (phase - gLastPhase);
    gLastPhase = phase;
    
    tmp = tmp - (0:length(gFFTworksp)-1)*expct;
    
%     qpd = floor(tmp2/pi);
%     for i=1:length(qpd)
%         if qpd(i) >= 0 && mod(qpd(i),2)==0
%             qpd(i) = qpd(i)+1;
%         elseif qpd(i) < 0 && mod(qpd(i),2)==0
%                 qpd(i) = qpd(i)-1;
%         end
%     end
%     tmp2 = tmp2- pi*qpd;
    tmp2=angle(complex(cos(tmp), sin(tmp)));
    
    tmp3 = osamp*tmp2/(2*pi);
    
    tmp4 = (0:length(gFFTworksp)-1)*freqPerBin + tmp3*freqPerBin;
    
    gAnaMagn = magn;
    gAnaFreq = tmp4;
    
    
    gSynMagn = zeros(1,length(gAnaMagn));
    gSynFreq = zeros(1,length(gAnaMagn));
    
    for k =1:length(gAnaMagn)
        index = floor((k-1)*pitchShift)+1;
        if index <= length(gAnaMagn)
            gSynMagn(index) = gSynMagn(index) + gAnaMagn(k); 
            gSynFreq(index) = gAnaFreq(k) * pitchShift; 
        end 
    end
    
    
    magn = gSynMagn;
    tmp = gSynFreq;
    
    tmp = tmp-(0:length(gFFTworksp)-1)*freqPerBin;
    
    tmp = tmp/freqPerBin;
    
    tmp = 2*pi*tmp/osamp;
    
    tmp = tmp + (0:length(gFFTworksp)-1)*expct;
    
    gSumPhase = gSumPhase + tmp;
    phase = gSumPhase;
    
    gFFTworksp = magn.* complex(cos(phase), sin(phase));
    
    gFFTworksp_real = ifft(gFFTworksp,fftFrameSize,'symmetric');
    
    windowed_frame = gFFTworksp_real.*window/(osamp*3/8);
    
    gOutputAccum = gOutputAccum+windowed_frame;
    
    outdata = gOutputAccum(1:stepSize);
    
    gOutputAccum = [gOutputAccum(stepSize+1:end), zeros(1,stepSize)];
    
    gInFIFO = [gInFIFO(stepSize+1:end), zeros(1,stepSize)];
    
    lastShift = pitchShift;
end