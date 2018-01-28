clear
vocpath = 'tests\gendata\';
d = dir([vocpath '*HPF.mat']);
fn = {d.name}';

fs = 44100;
delay = 1024-64;
nfreq = 49;
shift = -200:50:200;
rng = round(0.5*fs):round(2.0*fs);
Bhpf = fir1(1024,80/(fs/2),'high');

fdbkampAhpf = zeros(nfreq*2,length(shift),3,2); % freq x shift x algo x normalization
realampAhpf = zeros(nfreq*2,1); % freq x shift x algo x normalization
fdbkampChpf = zeros(nfreq*2,length(shift),3,2); % freq x shift x algo x normalization
realampChpf = zeros(nfreq*2,1); % freq x shift x algo x normalization
for f=1:length(fn)-2
    disp(fn{f}); tic;
    load([vocpath fn{f}]);
    % low cut for C (flat)
    filC = filtfilt(Bhpf,1,realHPF);
    realampChpf(f) = 20*log10(sqrt(mean(power(filC(rng),2))));
    % A-weighted filter
    filA = filterA(realHPF,fs);
    realampAhpf(f) = 20*log10(sqrt(mean(power(filA(rng),2))));
    for s=1:length(shift)
        for a=1:3
            for n=1:2
                temp = fdbkmatHPF(:,s,a,n);
                % low cut for C (flat)
                filC = filtfilt(Bhpf,1,temp);
                fdbkampChpf(f,s,a,n) = 20*log10(sqrt(mean(power(filC(rng+delay),2))));
                % A-weighted filter
                filA = filterA(temp,fs);
                fdbkampAhpf(f,s,a,n) = 20*log10(sqrt(mean(power(filA(rng+delay),2))));
            end
        end
    end
    disp(toc);
end
save tests/ampcalcHPF3.mat fdbkampAhpf realampAhpf fdbkampChpf realampChpf fs shift

