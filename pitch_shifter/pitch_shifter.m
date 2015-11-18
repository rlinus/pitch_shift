function [x_rec,x_plb,voice_onset] = pitch_shifter(mode,params)
%pitch_shifter is a tool for psychoacoustic pitch shift experiments
%   This function records audio and playbacks a pitch shifted version of it
%   in real time. A fixed sample rate and internal buffer size is used:
%       sampleRate: 44100
%       bufferSize: 64
%   All time parameters will therefore be rounded on integer multiples of
%   bufferSize/sampleRate seconds.
%
%   pitch_shifter(1,params):
%       Starts a trial (recording and playback). params is a struct with
%       the following fields:
%
%       shifterId:      Determines the pitch shift algorithm.
%                           0 = smbPitchShift [default]
%                           1 = cpvPitchShift
%                           2 = rubberband
%       deviceId:       Determines the audio device. Use print_device_list
%                       to get a list of possible IDs. Only ASIO devices
%                       are supported. [default: 0]
%       pitch_factor:   The pitch factor. [default: 1]
%       play_ref_sound: Play a reference sound at the beginning of the
%                       trial. [default: false]
%       ref_signal:     The reference signal. If not defined a synthesized
%                       e-piano sound will be played. [default: []]
%       ref_freq:       The frequency of the e-piano sound in Hz (only
%                       effective if ref_signal is empty). [default: 200]
%       ref_duration:   The duration of the e-piano sound in s (only
%                       effective if ref_signal is empty). [default: 1]
%       ref_amplitude:  The amplitude of the e-piano sound (only effective
%                       if ref_signal is empty). [default: 0.1]
%       voc_duration:   The duration of the trial in s, starting at the
%                       voice onset after the reference sound. [default: 1]
%       shift_onset:    The shift onset time in s relative to the voice
%                       onset time. [default: voc_duration/4]
%       shift_duration: The duration of the shift in s. [default: voc_duration/4]

end

