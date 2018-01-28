%PITCHLOOP plays a pitch shifted audio signal in a endless loop
%   PitchLoop(1, params, signal):
%       starts the playback of the pitch shifted signal in an endless loop.
%
%       signal:         the audio signal in vector form. The dynamic  
%                       range must be limited to +-1
%
%       params:    a struct with the following fields:
%
%       signal:         the audio signal in vector form. The dynamic  
%                       range must be limited to +-1
%       pitch_factor:   the pitch factor [default: 1]
%       shifterId:      Determines the pitch shift algorithm.
%                           0 = smbPitchShift [default]
%                           1 = cpvPitchShift
%                           2 = rubberband
%       windowSize:     The FFT window size of the Phase Vocoder.
%                           0: 512 [not available for rubberband]
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
%       rubberbandOptionPitch: pitch shift quality (only for rubberband).
%                       0:  OptionPitchHighSpeed - Use a method with a CPU
%                           cost that is relatively moderate and predictable.
%                           This may sound less clear than
%                           OptionPitchHighQuality, especially for large
%                           pitch shifts. 
%                       1:  OptionPitchHighQuality - Use the highest quality
%                           method for pitch shifting. This method has a
%                           CPU cost approximately proportional to the
%                           required frequency shift. [default]
%                       2:  OptionPitchHighConsistency - Use the method that
%                           gives greatest consistency when used to create
%                           small variations in pitch around the 1.0-ratio
%                           level. Unlike the previous two options, this
%                           avoids discontinuities when moving across the
%                           1.0 pitch scale in real-time; it also consumes
%                           more CPU than the others in the case where the
%                           pitch scale is exactly 1.0.
%       rubberbandOptionPhase: (only for rubberband)
%                       0:  OptionPhaseLaminar - Adjust phases when
%                           stretching in such a way as to try to retain the
%                           continuity of phase relationships between
%                           adjacent frequency bins whose phases are
%                           behaving in similar ways. [defualt]
%                       1:  OptionPhaseIndependent - Adjust the phase in
%                           each frequency bin independently from its
%                           neighbours. This usually results in a slightly
%                           softer, phasier sound.
%                       
%   PitchLoop(0, pitch_factor):
%       changes the pitch_factor
%
%   PitchLoop(-1):
%       stops the playback


