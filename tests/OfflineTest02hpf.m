clear
vocpath = '.\tests\vocoded\';
d = dir([vocpath '*.wav']);
fn = {d.name}';


fs = 44100;
shift = -200:50:200;
shift_fac = 2.^(shift./1200);
algo = [0 1 2]; % 0: smbPitchShift [default] / 1: cpvPitchShift / 2: rubberband
Bhpf = fir1(1024,80/(fs/2),'high');

% -------- shifter parameter setting -----------------------
frameSize = 64;
params.deviceId = 0;  
params.windowSize = 1; % 0: 512 / 1: 1024 [default] / 2: 2048 (fixed to 1024 if rubberband)
params.play_ref_sound = 0;
params.shift_full_trial = 1;
params.add_pink_noise = 0;
params.do_var = 0;
params.do_var_full_trial = 0;
params.feedback_gain = 1;
% -------------------------------------------------

for f=1:length(fn)
%for f=10:11
    disp(fn{f});tic;
    real = audioread([vocpath fn{f}]);
    wlen = length(real);
    wlen = floor(wlen/frameSize)*frameSize;
    real = real(1:wlen);
    params.voc_duration = wlen/fs;
    realHPF = filtfilt(Bhpf,1,real);
    fdbkmatHPF = zeros(wlen,length(shift),3,2); % shift x algo x normalization
    for n=1:2
        params.volume_normalization = n-1; 
        for a=1:3
            params.shifterId = algo(a); 
            for s=1:length(shift)
                params.pitch_factor = shift_fac(s);
                PsychPitchShifter(2,params,realHPF);
                [~,fdbk] = PsychPitchShifter(-1);
                fdbkmatHPF(:,s,a,n) = fdbk;
            end
        end
    end
    [a,b,c] = fileparts(fn{f});
    save(['tests\gendata\' b '_HPF.mat'],'fdbkmatHPF','realHPF');
    disp(toc);
end


