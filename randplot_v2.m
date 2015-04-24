%data = pitch_exp_recording;

Fs = data.Fs;

%% prepare and convert data for RTpitchMathworksRand
lev0 = find(data.pitch_levels_cents>0,1);
cents_pitches = [data.pitch_levels_cents(1:lev0-1) 0 data.pitch_levels_cents(lev0:end)];
h=1:length(cents_pitches);
levs = h(h~=lev0);


Levs = zeros(length(data.session_starts),1);
for i=1:length(Levs)
    Levs(i) = find(cents_pitches == data.pitch_level_sqs_cents(i),1);
end


cs = length(data.session_starts)-1;
Ys0 = cell(cs,1);
Ys = Ys0;
for i=1:cs
    Ys0{i} = data.y_r(data.session_starts(i):data.session_starts(i+1)-1);
    Ys{i} = data.y_ps(data.session_starts(i):data.session_starts(i+1)-1);
end
%Ys0{length(session_starts)} = y_r(session_starts(end):end);
%Ys{length(session_starts)} = y_ps(session_starts(end):end);


WindowLen = 256;
AnalysisLen = 64;

Vars=zeros(size(Levs));
do_var = 0;
computed_pitch=0;


RTpitchMathworksRandplot_log2
