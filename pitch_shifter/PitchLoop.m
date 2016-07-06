%PITCHLOOP plays a pitch shifted audio signal in a endless loop
%   PitchLoop(1, signal, pitch_factor, shifterId, deviceId, volume_normalization, tau, delta, windowSize):
%       starts the playback of the pitch shifted signal in an endless loop
%
%       signal:         the audio signal in vector form. The dynamic  
%                       range must be limited to +-1
%       pitch_factor:   the pitch factor [default: 1]
%       shifterId:      Determines the pitch shift algorithm.
%                           0 = smbPitchShift [default]
%                           1 = cpvPitchShift
%                           2 = rubberband
%       windowSize:     The FFT window size of the Phase Vocoder. For the
%                       rubberband algorithm the windowSize is not
%                       adjustable and fixed to 1024.
%                           0: 512
%                           1: 1024 [default]
%                           2: 2048
%       deviceId:       Determines the output audio device. Use
%                       PrintDeviceList to get a list of possible IDs.
%                       Only ASIO devices are supported. [default: 0]
%       volume_normalization: Normalizes the output signal, by applying the
%                       following transformation to the output signal s[t]:
%                       x[t] = (1-tau)* x[t-1]+tau*s[t]^2
%                       o[t] = s[t]*ln(sqrt(x[t])/delta+1)/(sqrt(x[t])/delta)
%                       [default: 0]
%       tau:            Volume normalisation constant (0<=tau<=1) [default: 0.5]
%       delta:          Volume normalisation constant (delta > 0) [default: 0.5]
%                       
%   PitchLoop(0, pitch_factor):
%       changes the pitch_factor
%
%   PitchLoop(-1):
%       stops the playback


