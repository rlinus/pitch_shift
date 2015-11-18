function pitch_loop(mode, signal, pitch_factor, shifterId, deviceId)
%pitch_loop plays pitch shifted audio signal in endless loop
%   pitch_loop(1, signal, pitch_factor, shifterId, deviceId):
%       starts the playback of the pitch shifted signal in an endless loop
%
%       signal:         the audio signal in vector form. The dynamic  
%                       range must be limited to +-1
%       pitch_factor:   the pitch factor [default: 1]
%       shifterId:      Determines the pitch shift algorithm.
%                           0 = smbPitchShift [default]
%                           1 = cpvPitchShift
%                           2 = rubberband
%       deviceId:       Determines the output audio device. Use
%                       print_device_list to get a list of possible IDs.
%                       Only ASIO devices are supported. [default: 0]
%                       
%   pitch_loop(0, pitch_factor):
%       changes the pitch_factor
%
%   pitch_loop(-1):
%       stops the playback

end

