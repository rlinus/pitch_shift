function outdata = perc_adapt_plot(data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for i=1:data.n_trials
    b = round((data.voice_onset_s(i)+0.3)*data.Fs/data.frameSize);
    e = round(0.1*data.Fs/data.frameSize);
    p = data.detected_pitch{i}(b:end-e);
    outdata.produced_freq_mean(i) = mean(p(p~=0));
end

x = 1:75;

%all trials
figure; hold on;
plot(x,data.pitch_level_sqs_cents);
plot(x,data.pitch_level_var_sqs_cents);
plot(x,data.perceived_produced_pitch);
plot(x,data.kind_of_trials*-100);
legend('mean applied shift','var applied shift', 'PPD')
title('all trials');

%feedback trials
figure; hold on;
x1 = 1:sum(data.kind_of_trials==1);
plot(x1,data.pitch_level_sqs_cents(data.kind_of_trials==1));
plot(x1,data.pitch_level_var_sqs_cents(data.kind_of_trials==1));
plot(x1,data.perceived_produced_pitch(data.kind_of_trials==1));
legend('mean applied shift','var applied shift', 'PPD')
title('feedback trials');

%motor catch trial (absolute performance)
figure; hold on;
x1 = 1:sum(data.kind_of_trials==2);
plot(x1,1200*log2(data.ref_freq_sqs(data.kind_of_trials==2)/data.ref_freq));
plot(x1,1200*log2(outdata.produced_freq_mean(data.kind_of_trials==2)/data.ref_freq));
plot(x1,data.pitch_level_sqs_cents(data.kind_of_trials==2));
title('motor catch trials: absolute performance');


%motor catch trial
x1 = 1:sum(data.kind_of_trials==2);
x2=x(data.kind_of_trials==2);
prevs = nan(1,length(x2));
for i=1:length(x2)
    for j=x2(i)-1:-1:1
        if data.kind_of_trials(j)==1
            prevs(i)=j;
            break;
        end
    end
end
x3=x1(~isnan(prevs));
prevs=prevs(~isnan(prevs));


figure; hold on;

plot(x1,1200*log2(outdata.produced_freq_mean(data.kind_of_trials==2)./data.ref_freq_sqs(data.kind_of_trials==2)));
plot(x1,data.pitch_level_sqs_cents(data.kind_of_trials==2));
plot(x1,data.perceived_produced_pitch(data.kind_of_trials==2));
plot(x3,data.pitch_level_var_sqs_cents(prevs));
title('motor catch trials');

%perceptual probe trials
x1 = 1:sum(data.kind_of_trials==3);
x2=x(data.kind_of_trials==3);
prevs = nan(1,length(x2));
for i=1:length(x2)
    for j=x2(i)-1:-1:1
        if data.kind_of_trials(j)==1
            prevs(i)=j;
            break;
        end
    end
end
x3=x1(~isnan(prevs));
prevs=prevs(~isnan(prevs));


figure; hold on;

plot(x1,data.pitch_level_sqs_cents(data.kind_of_trials==3));
plot(x1,data.perceived_produced_pitch(data.kind_of_trials==3));
plot(x3,data.pitch_level_var_sqs_cents(prevs));
title('perceptual probe trials');

end

