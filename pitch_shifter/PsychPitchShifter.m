%PSYCHPITCHSHIFTER is a psychophysiological pitch shift experiment tool
%   This function records audio and playbacks a pitch shifted version of it
%   in real time. A fixed sample rate and internal audio buffer size is
%   used:
%       sampleRate: 44100
%       bufferSize: 64
%   All time parameters will therefore be rounded on integer multiples of
%   bufferSize/sampleRate seconds.
%
%   PsychPitchShifter(1,params):
%       Starts a trial (recording and playback). params is a struct with
%       the following fields:
%
%       shifterId:      Determines the pitch shift algorithm.
%                           0: smbPitchShift [default]
%                           1: cpvPitchShift
%                           2: rubberband
%       windowSize:     The FFT window size of the Phase Vocoder.
%                           0: 512 [not available for rubberband]
%                           1: 1024 [default]
%                           2: 2048
%       deviceId:       Determines the audio device. Use the
%                       PrintDeviceList auxiliary function to get a list of
%                       possible IDs. Only ASIO devices are supported.
%                       [default: 0]
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
%       shift_duration: The duration of the shift in s.
%                       [default: voc_duration/4]
%       shift_full_trial: The shift is applied during the full trial.
%                       [default: 0]
%       do_var:         Makes the pitch shift variable, it follows a
%                       realization of a low-pass filtered white noise 
%                       stochastic process. [default: 0]
%       do_var_full_trial: Makes the pitch variable during the whole trial,
%                       not only during the shift. [default: 0]
%       std_dev:        the standard deviation of the white noise in cents.
%                       [defualt: 100]
%       fc:             the normalized cutoff frequency of the noise 
%                       filter. Can have values between 0 and 1. Set it to
%                       1 to bybass the filter. [defualt: 0.01]
%       T_var:          Sets the sampling time of the stochastic process as
%                       a multiple of the frame rate
%                       (=sampleRate/frameSize). [default: 1]
%       var_quant_size: Quantizies the amplitude of the shifts to prevent
%                       very small changes in pitch. Set it to 0 to turn
%                       this off. [default: 0]
%       do_control:     Uses a PI controller that outputs an additional
%                       shift as an apttempt to drive the produced pitch
%                       to full compensation. [default: 0]
%       kp:             kp Parameter of the PI controller. [default: 0]
%       ki:             ki Parameter of the PI controller. [default: 1]
%       control_ref_freq: The reference frequency for th PI controller
%                       [default: ref_freq]
%       add_pink_noise: Adds Pink Noise to the output signal. [default: 0]
%       noise_gain:     The amplitude of the noise. [default: 0.005]
%       adaptive_noise_lvl: makes the noise aplitude adaptive to the
%                       input sound volume with the gain determined by 
%                       noise_gain. [default: false]
%       min_noise_lvl:  the minimum noise level in the adaptive noise mode.
%                       [default: 0.005]
%       max_noise_lvl:  the maximum noise level in the adaptive noise mode.
%                       [default: 0.05]
%       feedback_gain:  The amplitude of the output voice signal. 
%                       [default: 1]
%       start_threshold: The threshold for the voice onset detection.
%                       [default: 0.01]
%       stop_threshold: The threshold for the voice offset detection.
%                       [default: 0.01]
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
%   state = PsychPitchShifter(0):
%       Returns a value according to the current state.
%           0:  the function returns 0 till voc_duration seconds after the
%               voice onset.
%           1:  the function returns 1 after the end of the trial (which is
%               voc_duration seconds after the voice onset), if the
%               voice amplitude is still greater than stop_threshold.
%           2:  the function returns 2 after the end of the trial, if the
%               voice amplitude has reached a value smaller than
%               stop_threshold.
%
%   [x,y,voice_on,stc_pf,var_pf,ctrl_pf,dpitch] = PsychPitchShifter(-1):
%       Stops the audio recording and playback and returns the data.
%           x:          The recorded audio data.
%           y:          The pitch shifted audio data.
%           voice_on:   The voice onset time in seconds.
%           stc_pf:     The sequence (one sample per frame) of static
%                       pitch factors.
%           var_pf:     The sequence of pitch factor fractions produced by
%                       the do_var feature.
%           ctrl_pf:    The sequence of pitch factor fractions produced by
%                       the do_control feature. Multiply all 3 factors to
%                       get the accumulated applied pitch shift.
%           dpitch:     The sequence of detected pitches. The do_control
%                       feature relies on pitch detection, this is the
%                       output of the internal pitch detection algorithm.
%
%    PsychPitchShifter(2,params,input):
%       Offline mode. Does not record from audio device but process a given
%       input signal. The call blocks till the whole signal is processed.
%       The processed signal can then be retrieved with:
%       [x,y,voice_on,stc_pf,var_pf,ctrl_pf,dpitch] = PsychPitchShifter(-1);
%       
%           input:      Audio signal as a double vector.
%           params:     Parameters as in PsychPitchShifter(1,params) (see
%                       above)
%__________________________________________________________________________
% (c) 2015 Linus Ruettimann (linus[dot]ruettimann[at]gmail[dot]com)
