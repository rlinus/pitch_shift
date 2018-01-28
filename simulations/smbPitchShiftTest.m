sampleRate = 44100;
frameSize = 64;
windowSize = 1024;

pitchShift=2^(0/1200);

T=5;
T_f = floor(T*sampleRate/frameSize);
T= T_f*frameSize/sampleRate;

t=(1/sampleRate:1/sampleRate:T);

f=320.2;
sig = sin(2*pi*f*t);%+0.5*sin(2*pi*2*f*t)+0.25*sin(2*pi*4*f*t);
%sig = data.y_r{2}(5000:end);

outdata=zeros(1,length(sig));
windowed_frame = zeros(T_f,windowSize);
for i=1:T_f
%     cents = -1*floor(i/100);
%     pitchShift=2^(cents/1200);
%pitchShift = 0.9+floor(i/100)*0.01;
    if i==200,pitchShift=2^(-200/1200);end;
    if i==1000,pitchShift=2^(-000/1200);end;
    if i==1500,pitchShift=2^(000/1200);end;
    %[outdata(1+(i-1)*frameSize:i*frameSize),windowed_frame(i,:)] = smbPitchShift(pitchShift,sig(1+(i-1)*frameSize:i*frameSize),i==1);
    [outdata(1+(i-1)*frameSize:i*frameSize),windowed_frame(i,:)] = cpvPitchShift(pitchShift,sig(1+(i-1)*frameSize:i*frameSize),i==1);
end

t=1:length(sig);
figure(1);
plot(t,sig,t,outdata);
figure(3);
% 
te = (1:windowSize);
k = 1100;
plot(te,windowed_frame(k,:),te+frameSize,windowed_frame(k+1,:),te+2*frameSize,windowed_frame(k+2,:),te+3*frameSize,windowed_frame(k+3,:));
% 
% b = windowed_frame(2100,:);
% plot(te,a,te,b);
%plot(te,windowed_frame(100,:),te,windowed_frame(1100,:));