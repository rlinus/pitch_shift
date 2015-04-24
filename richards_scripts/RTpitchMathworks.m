ratio=1.5;
WindowLen = 256*ratio;
AnalysisLen = .75*64*ratio;
SynthesisLen = round(65*ratio);
WinL=round(WindowLen/SynthesisLen*AnalysisLen);
Hopratio = SynthesisLen/AnalysisLen;

cent_pitch=1200*log(SynthesisLen/AnalysisLen);

windowSize=WindowLen; noverlap=AnalysisLen;

offline=1;
if offline
hAudioSource = dsp.AudioFileReader(...
  which('speech_dft_8kHz.wav'), ...
  'SamplesPerFrame', AnalysisLen, ...
  'OutputDataType', 'double');
else
hAudioSource = dsp.AudioRecorder('NumChannels',1,'SampleRate',Fs,'OutputDataType','double','SamplesPerFrame',noverlap,'BufferSizeSource','Property','BufferSize',512,'QueueDuration',1024/44100);
end

hbuf = dsp.Buffer(WindowLen, WindowLen - AnalysisLen);

% Create a Window System object, which is used for the ST-FFT. This object
% applies a window to the buffered input data.
hwin = dsp.Window('Hanning', 'Sampling', 'Periodic');

hfft = dsp.FFT;

% Create an IFFT System object, which is used for the IST-FFT.
hifft = dsp.IFFT('ConjugateSymmetricInput', true, ...
  'Normalize', false);


Fs = 8000;
%Fs=44100;
if offline
    hAudioOut = dsp.AudioPlayer('SampleRate', Fs);
else
    hAudioOut = dsp.AudioPlayer('SampleRate',Fs,'BufferSizeSource','Property','BufferSize',512,'QueueDuration',1024/Fs);
end
% Create a System object to log your data.
hslg = dsp.SignalSink;


%yprevwin = zeros(WindowLen-SynthesisLen, 1);
yprevwin = zeros(WinL-AnalysisLen, 1);
gain = 1/(WindowLen*sum(hanning(WindowLen,'periodic').^2)/SynthesisLen);
unwrapdata = 2*pi*AnalysisLen*(0:WindowLen-1)'/WindowLen;
yangle = zeros(WindowLen, 1);
firsttime = true;
%ly=0;
%yy=[];
%while ~isDone(hAudioSource)
yy=yangle(1:WindowLen/2);ylim([0 12]);
figure(1);clf; hh=semilogy(yy);drawnow
for ii=1:500
    y = step(hAudioSource);

 %   step(hAudioOut, y);    % Play back original audio

    % ST-FFT
    % FFT of a windowed buffered signal
    yfft = step(hfft, step(hwin, step(hbuf, y)));

    % Convert complex FFT data to magnitude and phase.
    ymag       = abs(yfft);
    yprevangle = yangle;
    yangle     = angle(yfft);
    yy=.9*yy+.1*ymag(1:WindowLen/2);
 %   ymag(400:500)=1e-3;
    if rem(ii,10)==0
%set(hh,'YData',yy); drawnow
    end
    % Synthesis Phase Calculation
    % The synthesis phase is calculated by computing the phase increments
    % between successive frequency transforms, unwrapping them, and scaling
    % them by the ratio between the analysis and synthesis hop sizes.
    yunwrap = (yangle - yprevangle) - unwrapdata;
    yunwrap = yunwrap - round(yunwrap/(2*pi))*2*pi;
    yunwrap = (yunwrap + unwrapdata) * Hopratio;
    if firsttime
        ysangle = yangle;
        firsttime = false;
    else
        ysangle = ysangle + yunwrap;
    end

    % Convert magnitude and phase to complex numbers.
    ys = ymag .* complex(cos(ysangle), sin(ysangle));

    % IST-FFT
    ywin  = step(hwin, step(hifft,ys));    % Windowed IFFT
    ywin=interpft(ywin,WinL);
    % Overlap-add operation
    figure(1);clf;plot(ywin(1:end-AnalysisLen,:),'b'); hold on;plot(yprevwin,'r');
    pause
    olapadd  = [ywin(1:end-AnalysisLen,:) + yprevwin; ...
        ywin(end-AnalysisLen+1:end,:)];
    yistfft  = olapadd(1:AnalysisLen,:);
    yprevwin = olapadd(AnalysisLen+1:end,:);

%     % Overlap-add operation
%     olapadd  = [ywin(1:end-SynthesisLen,:) + yprevwin; ...
%                 ywin(end-SynthesisLen+1:end,:)];
%     yistfft  = olapadd(1:SynthesisLen,:);
%     yprevwin = olapadd(SynthesisLen+1:end,:);

    % Compensate for the scaling that was introduced by the overlap-add
    % operation
    yistfft = yistfft * gain;
    step(hAudioOut,[yistfft yistfft]);
 %   step(hslg, yistfft);     % Log signal
end


release(hAudioSource);
release(hAudioOut);

break

loggedSpeech = hslg.Buffer(200:end)';
% Play time-stretched signal
p2 = audioplayer(loggedSpeech, Fs);
disp('Playing time-stretched signal...');
playblocking(p2);

% Play pitch-scaled signal
% The pitch-scaled signal is the time-stretched signal played at a higher
% sampling rate which produces a signal with a higher pitch.
Fs_new = Fs*(SynthesisLen/AnalysisLen);
p3 = audioplayer(loggedSpeech, Fs_new);
disp('Playing pitch-scaled signal...');
playblocking(p3);
