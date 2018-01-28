% show data
clear
load tests/ampcalcHPF

A2 = 110;
cent = 300:50:2700;
freq = A2*2.^(cent./1200);
adAM = fdbkampAhpf(1:49,:,:,:) - repmat(realampAhpf(1:49),[1 9 3 2]);
adAF = fdbkampAhpf(50:98,:,:,:) - repmat(realampAhpf(50:98),[1 9 3 2]);

yl = [-18 2];

fig = figure;
set(fig,'name','malevoice');
subplot(3,3,1)
plot(shift,adAM(:,:,1,1)','k',shift,mean(adAM(:,:,1,1)),'g',[-200 200],[0 0],'k:')
title('algorithm=1 / normalize=0')
set(gca,'tickdir','out','ylim',yl,'xtick',shift);
xlabel('Shift (cent)'); ylabel({'Loudness change (dB)','A-weighted'});

subplot(3,3,2)
plot(shift,adAM(:,:,2,1)','k',shift,mean(adAM(:,:,2,1)),'g',[-200 200],[0 0],'k:')
title('algorithm=2 / normalize=0')
set(gca,'tickdir','out','ylim',yl,'xtick',shift);
xlabel('Shift (cent)'); ylabel({'Loudness change (dB)','A-weighted'});

subplot(3,3,3)
plot(shift,adAM(:,:,3,1)','k',shift,mean(adAM(:,:,3,1)),'g',[-200 200],[0 0],'k:')
title('algorithm=3 / normalize=0')
set(gca,'tickdir','out','ylim',yl,'xtick',shift);
xlabel('Shift (cent)'); ylabel({'Loudness change (dB)','A-weighted'});

subplot(3,3,4)
plot(shift,adAM(:,:,1,2)','k',shift,mean(adAM(:,:,1,2)),'g',[-200 200],[0 0],'k:')
title('algorithm=1 / normalize=1')
set(gca,'tickdir','out','ylim',yl,'xtick',shift);
xlabel('Shift (cent)'); ylabel({'Loudness change (dB)','A-weighted'});

subplot(3,3,5)
plot(shift,adAM(:,:,2,2)','k',shift,mean(adAM(:,:,2,2)),'g',[-200 200],[0 0],'k:')
title('algorithm=2 / normalize=1')
set(gca,'tickdir','out','ylim',yl,'xtick',shift);
xlabel('Shift (cent)'); ylabel({'Loudness change (dB)','A-weighted'});

subplot(3,3,6)
plot(shift,adAM(:,:,3,2)','k',shift,mean(adAM(:,:,3,2)),'g',[-200 200],[0 0],'k:')
title('algorithm=3 / normalize=1')
set(gca,'tickdir','out','ylim',yl,'xtick',shift);
xlabel('Shift (cent)'); ylabel({'Loudness change (dB)','A-weighted'});

subplot(3,3,7);
plot(freq,adAM(:,[2 4 6 8],1,2),'-',[0 600],[0 0],'k:')
title('algorithm=1 / normalize=1')
legend({'-150c','-50c','+50c','+150c'},'location','best');
set(gca,'tickdir','out','ylim',yl,'xlim',[100 550]);
xlabel('Original Fo (Hz)'); ylabel({'Loudness change (dB)','A-weighted'});

subplot(3,3,8);
plot(freq,adAM(:,[2 4 6 8],2,2),'-',[0 600],[0 0],'k:')
title('algorithm=2 / normalize=1')
legend({'-150c','-50c','+50c','+150c'},'location','best');
set(gca,'tickdir','out','ylim',yl,'xlim',[100 550]);
xlabel('Original Fo (Hz)'); ylabel({'Loudness change (dB)','A-weighted'});

subplot(3,3,9);
plot(freq,adAM(:,[2 4 6 8],3,2),'-',[0 600],[0 0],'k:')
title('algorithm=3 / normalize=1')
legend({'-150c','-50c','+50c','+150c'},'location','best');
set(gca,'tickdir','out','ylim',yl,'xlim',[100 550]);
xlabel('Original Fo (Hz)'); ylabel({'Loudness change (dB)','A-weighted'});

%%
fig = figure;
set(fig,'name','femalevoice');

subplot(3,3,1)
plot(shift,adAF(:,:,1,1)','k',shift,mean(adAF(:,:,1,1)),'g',[-200 200],[0 0],'k:')
title('algorithm=1 / normalize=0')
set(gca,'tickdir','out','ylim',yl,'xtick',shift);
xlabel('Shift (cent)'); ylabel({'Loudness change (dB)','A-weighted'});

subplot(3,3,2)
plot(shift,adAF(:,:,2,1)','k',shift,mean(adAF(:,:,2,1)),'g',[-200 200],[0 0],'k:')
title('algorithm=2 / normalize=0')
set(gca,'tickdir','out','ylim',yl,'xtick',shift);
xlabel('Shift (cent)'); ylabel({'Loudness change (dB)','A-weighted'});

subplot(3,3,3)
plot(shift,adAF(:,:,3,1)','k',shift,mean(adAF(:,:,3,1)),'g',[-200 200],[0 0],'k:')
title('algorithm=3 / normalize=0')
set(gca,'tickdir','out','ylim',yl,'xtick',shift);
xlabel('Shift (cent)'); ylabel({'Loudness change (dB)','A-weighted'});

subplot(3,3,4)
plot(shift,adAF(:,:,1,2)','k',shift,mean(adAF(:,:,1,2)),'g',[-200 200],[0 0],'k:')
title('algorithm=1 / normalize=1')
set(gca,'tickdir','out','ylim',yl,'xtick',shift);
xlabel('Shift (cent)'); ylabel({'Loudness change (dB)','A-weighted'});

subplot(3,3,5)
plot(shift,adAF(:,:,2,2)','k',shift,mean(adAF(:,:,2,2)),'g',[-200 200],[0 0],'k:')
title('algorithm=2 / normalize=1')
set(gca,'tickdir','out','ylim',yl,'xtick',shift);
xlabel('Shift (cent)'); ylabel({'Loudness change (dB)','A-weighted'});

subplot(3,3,6)
plot(shift,adAF(:,:,3,2)','k',shift,mean(adAF(:,:,3,2)),'g',[-200 200],[0 0],'k:')
title('algorithm=3 / normalize=1')
set(gca,'tickdir','out','ylim',yl,'xtick',shift);
xlabel('Shift (cent)'); ylabel({'Loudness change (dB)','A-weighted'});

subplot(3,3,7);
plot(freq,adAF(:,[2 4 6 8],1,2),'-',[0 600],[0 0],'k:')
title('algorithm=1 / normalize=1')
legend({'-150c','-50c','+50c','+150c'},'location','best');
set(gca,'tickdir','out','ylim',yl,'xlim',[100 550]);
xlabel('Original Fo (Hz)'); ylabel({'Loudness change (dB)','A-weighted'});

subplot(3,3,8);
plot(freq,adAF(:,[2 4 6 8],2,2),'-',[0 600],[0 0],'k:')
title('algorithm=2 / normalize=1')
legend({'-150c','-50c','+50c','+150c'},'location','best');
set(gca,'tickdir','out','ylim',yl,'xlim',[100 550]);
xlabel('Original Fo (Hz)'); ylabel({'Loudness change (dB)','A-weighted'});

subplot(3,3,9);
plot(freq,adAF(:,[2 4 6 8],3,2),'-',[0 600],[0 0],'k:')
title('algorithm=3 / normalize=1')
legend({'-150c','-50c','+50c','+150c'},'location','best');
set(gca,'tickdir','out','ylim',yl,'xlim',[100 550]);
xlabel('Original Fo (Hz)'); ylabel({'Loudness change (dB)','A-weighted'});

%%
figure
%%
adswAM = fdbkampAhpf(99,:,:,:) - repmat(realampAhpf(99),[1 9 3 2]);
adswAF = fdbkampAhpf(100,:,:,:) - repmat(realampAhpf(100),[1 9 3 2]);
adswAM = squeeze(adswAM);
adswAF = squeeze(adswAF);
Fdat = load('sweep_female_F3-F4_HPF.mat');
Mdat = load('sweep_male_C3-C4_HPF.mat');

rfs = 1000;
win = kaiser(round(0.1*fs),20);
realswFenv = resample(10*log(filter(win,1,Fdat.realHPF.^2)),rfs,fs);
realswMenv = resample(10*log(filter(win,1,Mdat.realHPF.^2)),rfs,fs);
tempF = Fdat.fdbkmatHPF(:,:);
tempM = Mdat.fdbkmatHPF(:,:);
tempFr = resample(10*log(filter(win,1,tempF.^2)),rfs,fs);
tempMr = resample(10*log(filter(win,1,tempM.^2)),rfs,fs);
fdbkswFenv = reshape(tempFr,[3000 9 3 2]);
fdbkswMenv = reshape(tempMr,[3000 9 3 2]);
swdiffF = fdbkswFenv - repmat(realswFenv,[1 9 3 2]);
swdiffM = fdbkswMenv - repmat(realswMenv,[1 9 3 2]);

%%
rng = 520 + (1:2000);
swdMr = swdiffM(rng,:,:,:);
tvec = (1:length(rng))'/rfs;
jcm = jet;
co = jcm(2:7:64,:);
set(groot,'defaultAxesColorOrder',co)
yl = [-40 5];

fig = figure;
set(fig,'name','male sweep');
subplot(2,3,1)
plot(tvec*1000,swdMr(:,:,1,1));
set(gca,'tickdir','out','ylim',yl); 
title('algorithm=1 / normalize=0')
xlabel('Time (ms)'); ylabel({'Loudness difference (dB)','A-weighted'});

subplot(2,3,2)
plot(tvec*1000,swdMr(:,:,2,1));
set(gca,'tickdir','out','ylim',yl); 
title('algorithm=2 / normalize=0')
xlabel('Time (ms)'); ylabel({'Loudness difference (dB)','A-weighted'});

subplot(2,3,3)
plot(tvec*1000,swdMr(:,:,3,1));
set(gca,'tickdir','out','ylim',yl); 
title('algorithm=3 / normalize=0')
xlabel('Time (ms)'); ylabel({'Loudness difference (dB)','A-weighted'});

subplot(2,3,4)
plot(tvec*1000,swdMr(:,:,1,2));
set(gca,'tickdir','out','ylim',yl); 
title('algorithm=1 / normalize=1')
xlabel('Time (ms)'); ylabel({'Loudness difference (dB)','A-weighted'});

subplot(2,3,5)
plot(tvec*1000,swdMr(:,:,2,2));
set(gca,'tickdir','out','ylim',yl); 
title('algorithm=2 / normalize=1')
xlabel('Time (ms)'); ylabel({'Loudness difference (dB)','A-weighted'});

subplot(2,3,6)
plot(tvec*1000,swdMr(:,:,3,2));
set(gca,'tickdir','out','ylim',yl); 
title('algorithm=3 / normalize=1')
xlabel('Time (ms)'); ylabel({'Loudness difference (dB)','A-weighted'});
legend(cellstr(num2str(shift')))

%%
rng = 530 + (1:2000);
swdFr = swdiffF(rng,:,:,:);
fig = figure;
set(fig,'name','female sweep');
subplot(2,3,1)
plot(tvec*1000,swdFr(:,:,1,1));
set(gca,'tickdir','out','ylim',yl); 
title('algorithm=1 / normalize=0')
xlabel('Time (ms)'); ylabel({'Loudness difference (dB)','A-weighted'});

subplot(2,3,2)
plot(tvec*1000,swdFr(:,:,2,1));
set(gca,'tickdir','out','ylim',yl); 
title('algorithm=2 / normalize=0')
xlabel('Time (ms)'); ylabel({'Loudness difference (dB)','A-weighted'});

subplot(2,3,3)
plot(tvec*1000,swdFr(:,:,3,1));
set(gca,'tickdir','out','ylim',yl); 
title('algorithm=3 / normalize=0')
xlabel('Time (ms)'); ylabel({'Loudness difference (dB)','A-weighted'});

subplot(2,3,4)
plot(tvec*1000,swdFr(:,:,1,2));
set(gca,'tickdir','out','ylim',yl); 
title('algorithm=1 / normalize=1')
xlabel('Time (ms)'); ylabel({'Loudness difference (dB)','A-weighted'});

subplot(2,3,5)
plot(tvec*1000,swdFr(:,:,2,2));
set(gca,'tickdir','out','ylim',yl); 
title('algorithm=2 / normalize=1')
xlabel('Time (ms)'); ylabel({'Loudness difference (dB)','A-weighted'});

subplot(2,3,6)
plot(tvec*1000,swdFr(:,:,3,2));
set(gca,'tickdir','out','ylim',yl); 
title('algorithm=3 / normalize=1')
xlabel('Time (ms)'); ylabel({'Loudness difference (dB)','A-weighted'});
legend(cellstr(num2str(shift')))
