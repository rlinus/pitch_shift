function [noise] = get_pink_noise(num_samples)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
H = dsp.ColoredNoise('SamplesPerFrame', num_samples, 'OutputDataType', 'double');

noise = H.step;
end

