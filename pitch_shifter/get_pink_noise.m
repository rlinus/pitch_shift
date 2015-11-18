function [noise] = get_pink_noise(num_samples)
%this function is called by pitch_shifter and returns a vector of pink
%noise
H = dsp.ColoredNoise('SamplesPerFrame', num_samples, 'OutputDataType', 'double');

noise = H.step;
end

