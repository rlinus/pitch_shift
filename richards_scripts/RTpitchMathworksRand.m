%% Pitch shift ramdonly changes to up or down every even epoch (defined by pauses). On odd trials the pitch is zero

% Parameters
Amp_thresh=1; % above that sound amplitude no change in shift amount (parameter changes only during pauses)
Amp_thresh_high=5; % new idea: the true threshold is Amp_thresh_high, but Amp_thresh must be passed before

Min_Dur_per_level_s=1; % Pitch shift stays at least Min_Dur_per_level_s seconds at the specified level
T_s=20; % time in seconds

% lev0=2;
% levs=[1 3];

lev0=3;
levs=[1 2 4 5];
%levs=[1 3];

do_var=0;
std_lev=20;
varprob=1;

offline=0;


computed_pitch=0;

log_mult=16;
log_div=log_mult/4;

%% pitch sound hardware parameters
ratio=1;%1.5;
scal2=1;
WindowLen = 256*ratio;
AnalysisLen = scal2*64*ratio;
noverlap=AnalysisLen;
unwrapdata = 2*pi*AnalysisLen*(0:WindowLen-1)'/WindowLen;
windowSize=WindowLen;
yangle = zeros(WindowLen, 1);
firsttime = true;

SynthesisLens = round(scal2*round((64-2*lev0+2*[levs(levs<lev0) lev0 levs(levs>lev0)])*ratio));
cents_pitches=1200*log(SynthesisLens/AnalysisLen);



Fs=44100;
if offline
    hAudioSource = dsp.AudioFileReader(...
        which('blackriders_01_crane_64kb.wav'), ...
        'SamplesPerFrame', AnalysisLen, ...
        'OutputDataType', 'double');
else
    hAudioSource = dsp.AudioRecorder('NumChannels',1,'SampleRate',Fs,'OutputDataType','double','SamplesPerFrame',noverlap,'BufferSizeSource','Property','BufferSize',512,'QueueDuration',1024/Fs);
    %hAudioSource = dsp.AudioRecorder('DeviceName','MOTU Audio ASIO', 'NumChannels',1,'SampleRate',Fs,'OutputDataType','double','SamplesPerFrame',noverlap,'BufferSizeSource','Property','BufferSize',512,'QueueDuration',1024/Fs);
end
hbuf = dsp.Buffer(WindowLen, WindowLen - AnalysisLen);

% Create a Window System object, which is used for the ST-FFT. This object
% applies a window to the buffered input data.
hwin = dsp.Window('Hanning', 'Sampling', 'Periodic');
hfft = dsp.FFT;

% Create an IFFT System object, which is used for the IST-FFT.
hifft = dsp.IFFT('ConjugateSymmetricInput', true, ...
    'Normalize', false);
hAudioOut = dsp.AudioPlayer('SampleRate',Fs,'BufferSizeSource','Property','BufferSize',512,'QueueDuration',1024/Fs);
%hAudioOut = dsp.AudioPlayer('DeviceName','MOTU Audio ASIO', 'SampleRate',Fs,'BufferSizeSource','Property','BufferSize',512,'QueueDuration',1024/Fs);

% Create a System object to log your data.
hslg = dsp.SignalSink;

SynthesisLen = scal2*64*ratio; % first initialization = normal
WinL=round(WindowLen/SynthesisLen*AnalysisLen);
Hopratio = SynthesisLen/AnalysisLen;
cent_pitch=1200*log(SynthesisLen/AnalysisLen);

yprevwin = zeros(WinL-AnalysisLen, 1);
gain = 1/(WindowLen*sum(hanning(WindowLen,'periodic').^2)/SynthesisLen);
yy=yangle(1:WindowLen/2);ylim([0 12]);
figure(1);clf; hh=semilogy(yy);drawnow

T=ceil(T_s*Fs/AnalysisLen);
ymags=zeros(1,T);
yall=zeros(T*AnalysisLen,1);
yall2=yall;
Ys=cell(1,T_s*2); Ys0=Ys;
Ts=zeros(1,T_s*2);
Levs=Ts; Levs(1)=lev0;
Ts(1)=0;
Vars=zeros(size(Levs));
cs=0;    new_synt=0;
var_cumul=0; var_cumuls=zeros(1,T);
tic
delta=0.01;
new_synth=0;

% hh=zeros(1,T);
% for i=2:T
% hh(i)=(1-delta)*hh(i-1)+delta*std_lev*randn(1);
% end
% figure(1);clf;plot(hh)
%%
did_pass_low=0;
for ii=1:T
    y = step(hAudioSource);
    % y=y';
    %  y=y(:,1);
    yall((ii-1)*AnalysisLen+1:ii*AnalysisLen)=y;
    %   y=filter(DATfilt.aa,DATfilt.bb,y);
    
    %   step(hAudioOut, y);    % Play back original audio
    
    % ST-FFT
    % FFT of a windowed buffered signal
    yfft = step(hfft, step(hwin, step(hbuf, y)));
    
    % Convert complex FFT data to magnitude and phase.
    ymag       = abs(yfft);
    ymag(1:2)=0; ymag(end-1:end)=0;  ymag(WindowLen/2-WindowLen/8:WindowLen/2+WindowLen/8)=0; % high-pass filter
    
    ymags(ii)=sum(ymag); % sound amplitude
    ymag=log_mult*log(1+ymag/log_div); % logarithmic sound compression
    yprevangle = yangle;
    yangle     = angle(yfft);
    
    
    %    yy=.9*yy+.1*ymag(1:WindowLen/2); % running average spectrum for plotting
    %   ymag(400:500)=1e-3;
    %   if rem(ii,10)==0
    %      set(hh,'YData',yy); drawnow
    %   end
    
    % switch shift level
    if ii>1 && ymags(ii-1)<Amp_thresh && (ii-1)*AnalysisLen-Ts(cs+1)>(1+(Levs(cs+1)==lev0))*Min_Dur_per_level_s*Fs
        did_pass_low=1;
    end
    
    %if ii>1 && ymags(ii-1)<Amp_thresh && ymags(ii)>Amp_thresh && (ii-1)*AnalysisLen-Ts(cs+1)>Min_Dur_per_level_s*Fs %
    % new
    if did_pass_low && ymags(ii)>Amp_thresh_high %
        cs=cs+1;
        Ts(cs+1)=(ii-1)*AnalysisLen;
        Ys{cs}=yall2(Ts(cs)+1:Ts(cs+1)); % shifted pitch
        Ys0{cs}=yall(Ts(cs)+1:Ts(cs+1));
        if rem(cs,2)==0
            Levs(cs+1)=lev0; % normal
            var_cumul=0;
        else
            which_lev_i=ceil(length(levs)*rand(1));
            Levs(cs+1)=levs(which_lev_i);
            
            if do_var
                Vars(cs+1)=rand(1)>1-varprob;
            end
        end
        new_synth=1;
        did_pass_low=0;
    end
    
    if (~do_var && new_synth) || do_var %(do_var && rem(ii,4)==0)
        new_synt=0;
        var_cumul=(1-delta)*var_cumul+delta*Vars(cs+1)*std_lev*randn(1);
        var_cumuls(ii)=var_cumul;
        SynthesisLen = round(scal2*round((64-2*lev0+2*Levs(cs+1))*ratio)+var_cumul);
        WinLold=WinL;
        WinL=round(WindowLen/SynthesisLen*AnalysisLen);
        Hopratio = SynthesisLen/AnalysisLen;
        cent_pitch=1200*log(SynthesisLen/AnalysisLen);
        if WinL-AnalysisLen>length(yprevwin)
            yprevwin=[yprevwin; zeros(WinL-AnalysisLen-length(yprevwin),1)];
        else
            yprevwin=yprevwin(1:WinL-AnalysisLen);
        end
        gain = 1/(WindowLen*sum(hanning(WindowLen,'periodic').^2)/SynthesisLen);
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
    ywin0  = step(hwin, step(hifft,ys));    % Windowed IFFT
    %   ywin0=filter(DATfilt.aa,DATfilt.bb,ywin0);
    
    ywin=interpft(ywin0,WinL);
    % Overlap-add operation
    olapadd  = [ywin(1:end-AnalysisLen,:) + yprevwin; ...
        ywin(end-AnalysisLen+1:end,:)];
    yistfft  = olapadd(1:AnalysisLen,:);
    yprevwin = olapadd(AnalysisLen+1:end,:);
    % Compensate for the scaling that was introduced by the overlap-add
    % operation
    yistfft = yistfft * gain;
    yall2((ii-1)*AnalysisLen+1:ii*AnalysisLen)=yistfft;
    
    step(hAudioOut,[yistfft yistfft]);
    %   step(hslg, yistfft);     % Log signal
end
toc

release(hAudioSource);
release(hAudioOut);
figure(7); clf;plot(ymags); hold on; plot([1 length(ymags)],Amp_thresh*[1 1],'r');
figure(8);clf; [y,x]=hist(ymags,100);stairs(x,y,'LineWidth',0.8); hold on;   plot(x,log_mult*log(1+x/log_div),'r');
ylabel('num (log)'); xlabel('ymagnitudes');
ylim([0 max(x)]);% logarithmic sound compression

figure(9);clf; plot(var_cumuls);
%RTpitchMathworksRandplot
