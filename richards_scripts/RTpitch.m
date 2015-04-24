windowSize=1024; noverlap=256;
n=windowSize; hop=noverlap;
scf = 1.0;
r=.8;
%har = dsp.AudioRecorder('DeviceName','ASIO','NumChannels',4,'SampleRate',44100,'OutputDataType','double','SamplesPerFrame',windowSize - noverlap,'BufferSizeSource','Property','BufferSize',512,'QueueDuration',1024/44100);
% Handle to audio player
hap = dsp.AudioPlayer('SampleRate',44100,'BufferSizeSource','Property','BufferSize',512,'QueueDuration',1024/44100);
H = dsp.AudioRecorder('NumChannels',1,'SampleRate',44100,'OutputDataType','double','SamplesPerFrame',noverlap,'BufferSizeSource','Property','BufferSize',512,'QueueDuration',1024/44100);
N=10*44100/256;
a=zeros(noverlap,N);
aout=zeros(noverlap,N);
t=0:r:44100*10;
ct=0;
ncols=8;
ti=0;
for i=1:N
    audio = step(H); % Input from microphone
    a(:,i)=audio;
    ct=ct+1;
    if i>ncols+3
        x=a(:,i-ncols-3:i); x=x(:);
        %  x=resample(x,44100,round(44100*r));
        %  x=interpft(x,round(length(x)*r));
        % x=x(end-hop*(ncols):end);
        X = scf * stft(x, n, n, hop);
        [rows, cols] = size(X);
        ti_old=ti;
        ti=t(t>ct-ncols+4 & t<ct);% t = 0:r:(cols-2);
        %         if ti(end)
        %         end
        X2 = pvsample(X, ti-(ct-ncols+4), hop);
        y = istft(X2, n, n, hop)';
        %aout(:,i)=y(end-2*256+1:end-1*256);
        ao=y(end-4*hop+1:end-3*hop);
        aout(:,i)=ao;
        step(hap,100*[ao ao]);
    end
end
pause(H.QueueDuration);  % Wait until audio plays to the end

release(H); % close audio input device, release resources

release(hap);  % close audio output device, release resources
%figure(1);clf;plot(a(:))
%soundsc(aout(:),44100/r);

break
n=256*4; hop=256;
X = scf * stft(a(:), n, n, hop);

[rows, cols] = size(X);
t = 0:r:(cols-2);
% Have to stay two cols off end because (a) counting from zero, and
% (b) need col n AND col n+1 to interpolate

% Generate the new spectrogram
X2 = pvsample(X, t, hop);

% Invert to a waveform
y = istft(X2, n, n, hop)';