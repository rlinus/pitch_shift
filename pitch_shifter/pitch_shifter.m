%PITCH_SHIFTER is a tool for psychoacoustic pitch shift experiments
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
%       shift_full_trial: The shift is applied during the full trial.
%                       [default: 0]
%       do_var:         Makes the pitch shift variable, it follows a
%                       realization of a low-pass filtered white noise 
%                       stochastic process. [default: 0]
%       do_var_full_trial: Makes the pitch variable during the whole trial,
%                       not only during the shift. [default: 0]
%       std_dev:        the standard deviation of the white noise in cents.
%                       [defualt: 100]
%       fc:             the cutoff frequency of the noise filter normalized
%                       to 1. Set it to 1 to bybass the filter. [defualt: 0.01]
%       T_var:          Sets the sampling time of the stochastic process as
%                       a multiple of the frame rate
%                       (=sampleRate/frameSize). [default: 1]
%       var_quant_size: Quantizies the amplitude of the shifts to prevent
%                       very small changes in pitch. Set it to 0 to turn
%                       this off. [default: 0]
%       do_control:     Uses a PI controller that outputs an additional
%                       shift as an apttempt to drive the produced pitch
%                       to full compensation. [default: 0]
%       kp:             Parameter of the PI controller. [default: 0]
%       ki:             Parameter of the PI controller. [default: 1]
%       add_pink_noise: Adds Pink Noise to the output signal. [default: 0]
%       noise_gain:     The amplitude of the noise. [default: 0.005]
%       feedback_gain:  The amplitude of the output voice signal. 
%                       [default: 1]
%       start_threshold: The threshold for the voice onset detection.
%                       [default: 0.01]
%       stop_threshold: The threshold for the voice offset detection.
%                       [default: 0.01]
%
%   state = pitch_shifter(0):
%       Returns a value according to the current state.
%           0:  the function returns 0 till voc_duration second after the
%               voice onset.
%           1:  the function returns 1 after the end of the trial, if the
%               voice amplitude is greater than stop_threshold.
%           2:  the function returns 2 after the end of the trial, if the
%               voice amplitude has reached a value smaller than
%               stop_threshold.
%
%   [x,y,voice_onset,static_pitch_factor, var_pitch_factor, control_pitch_factor, detected_pitch] = pitch_shifter(-1):
%       Stops the audio recording and playback and returns the data.
%           x:              The recorded audio data.
%           y:              The pitch shifted audio data.
%           voice_onset:    The voice onset time in seconds.
%           static_pitch_factor:    The sequence (one sample per frame)
%                                   of static pitch factors.
%           var_pitch_factor:       The sequence of pitch factor fractions
%                                   produced by the do_var feature.
%           control_pitch_factor:   The sequence of pitch factor fractions
%                                   produced by the do_control feature.
%                                   Multiply all 3 factors to get the
%                                   applied pitch shift.
%           detected_pitch:         The sequence of detected pitches. The
%                                   do_control feature relies on pitch
%                                   detection, this is the output of the
%                                   internal pitch detection algorithm.


