function make(varargin)
    %RAWWAVEPATH = 'C:/Users/Linus/Documents/MATLAB/pitch_shift/rt_pitch_shifter/stk-4.5.0/rawwaves/';
    RAWWAVEPATH = 'E:/Data/Linus/GitHub/pitch_shift/rt_pitch_shifter/stk-4.5.0/rawwaves/';
    
    currentFolder = pwd;
    cd(fileparts(which(mfilename)));
    cd('bin');
    
    if isempty(varargin{1}) || ~ischar(varargin{1})
        fprintf('specify target');
        return;
    end
    
    path_param = ['-DRAWWAVEPATH=\"' RAWWAVEPATH '\"'];
    
    if strcmpi(varargin{1},'libs')
        % compile extern code
        mex ../stk-4.5.0/src/include/asio.cpp -DWIN32 -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include
        mex ../stk-4.5.0/src/include/asiolist.cpp  -DWIN32 -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include
        mex ../stk-4.5.0/src/include/asiodrivers.cpp  -DWIN32 -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include
        mex ../stk-4.5.0/src/include/iasiothiscallresolver.cpp  -DWIN32 -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include
        
        mex ../stk-4.5.0/src/Stk.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/RtAudio.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__  -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/FM.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/ADSR.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/FileLoop.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/SineWave.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/Wurley.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/FileWvIn.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/FileRead.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/TwoZero.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/Rhodey.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/BeeThree.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/Drummer.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/OnePole.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/LentPitShift.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/PitShift.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/Delay.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        mex ../stk-4.5.0/src/DelayL.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -I../stk-4.5.0/src/include -I../stk-4.5.0/include
        
        mex ../rubberband-1.8.1/src/RubberBandStretcher.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -I../rubberband-1.8.1

        mex ../rubberband-1.8.1/src/StretchCalculator.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -I../rubberband-1.8.1
        mex ../rubberband-1.8.1/src/StretcherChannelData.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -I../rubberband-1.8.1
        mex ../rubberband-1.8.1/src/StretcherImpl.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -I../rubberband-1.8.1
        mex ../rubberband-1.8.1/src/StretcherProcess.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -I../rubberband-1.8.1

        mex ../rubberband-1.8.1/src/audiocurves/SpectralDifferenceAudioCurve.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -I../rubberband-1.8.1 -I../rubberband-1.8.1/src
        mex ../rubberband-1.8.1/src/audiocurves/SilentAudioCurve.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -I../rubberband-1.8.1/src
        mex ../rubberband-1.8.1/src/audiocurves/PercussiveAudioCurve.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -I../rubberband-1.8.1/src
        mex ../rubberband-1.8.1/src/audiocurves/HighFrequencyAudioCurve.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -I../rubberband-1.8.1/src
        mex ../rubberband-1.8.1/src/audiocurves/ConstantAudioCurve.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -I../rubberband-1.8.1/src
        mex ../rubberband-1.8.1/src/audiocurves/CompoundAudioCurve.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -I../rubberband-1.8.1/src

        mex ../rubberband-1.8.1/src/base/Profiler.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -I../rubberband-1.8.1/src

        mex ../rubberband-1.8.1/src/dsp/Resampler.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -I../rubberband-1.8.1/src
        mex ../rubberband-1.8.1/src/dsp/FFT.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -I../rubberband-1.8.1/src
        mex ../rubberband-1.8.1/src/dsp/AudioCurveCalculator.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -I../rubberband-1.8.1 -I../rubberband-1.8.1/src

        mex ../rubberband-1.8.1/src/getopt/getopt_long.c -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -I../rubberband-1.8.1 -I../rubberband-1.8.1/src
        mex ../rubberband-1.8.1/src/getopt/getopt.c -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -I../rubberband-1.8.1 -I../rubberband-1.8.1/src

        mex ../rubberband-1.8.1/src/kissfft/kiss_fftr.c -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -I../rubberband-1.8.1 -I../rubberband-1.8.1/src
        mex ../rubberband-1.8.1/src/kissfft/kiss_fft.c -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -I../rubberband-1.8.1 -I../rubberband-1.8.1/src

        mex ../rubberband-1.8.1/src/speex/resample.c -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -I../rubberband-1.8.1 -I../rubberband-1.8.1/src

        mex ../rubberband-1.8.1/src/system/VectorOpsComplex.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -I../rubberband-1.8.1 -I../rubberband-1.8.1/src
        mex ../rubberband-1.8.1/src/system/Thread.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -I../rubberband-1.8.1 -I../rubberband-1.8.1/src
        mex ../rubberband-1.8.1/src/system/sysutils.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -I../rubberband-1.8.1 -I../rubberband-1.8.1/src
        mex ../rubberband-1.8.1/src/system/Allocators.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -I../rubberband-1.8.1 -I../rubberband-1.8.1/src
    elseif strcmpi(varargin{1},'rt_pitch_shifter')
        % compile custom code
        mex ../rt_pitch_shifter.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -g -c -I../stk-4.5.0/include -I../rubberband-1.8.1/rubberband

        % link everything
        mex rt_pitch_shifter.obj Stk.obj RtAudio.obj asio.obj asiolist.obj asiodrivers.obj iasiothiscallresolver.obj RubberBandStretcher.obj StretchCalculator.obj StretcherChannelData.obj StretcherImpl.obj StretcherProcess.obj SpectralDifferenceAudioCurve.obj SilentAudioCurve.obj PercussiveAudioCurve.obj HighFrequencyAudioCurve.obj ConstantAudioCurve.obj CompoundAudioCurve.obj Profiler.obj Resampler.obj FFT.obj AudioCurveCalculator.obj getopt_long.obj getopt.obj kiss_fftr.obj kiss_fft.obj resample.obj VectorOpsComplex.obj Thread.obj sysutils.obj Allocators.obj
    elseif strcmpi(varargin{1},'vowel_shifter_rubberband')
        mex ../vowel_shifter_rubberband.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -g -c -I../stk-4.5.0/include -I../rubberband-1.8.1/rubberband
        mex vowel_shifter.obj Stk.obj RtAudio.obj FM.obj ADSR.obj FileLoop.obj Delay.obj DelayL.obj SineWave.obj Wurley.obj Rhodey.obj BeeThree.obj Drummer.obj OnePole.obj LentPitShift.obj PitShift.obj FileWvIn.obj FileRead.obj TwoZero.obj asio.obj asiolist.obj asiodrivers.obj iasiothiscallresolver.obj RubberBandStretcher.obj StretchCalculator.obj StretcherChannelData.obj StretcherImpl.obj StretcherProcess.obj SpectralDifferenceAudioCurve.obj SilentAudioCurve.obj PercussiveAudioCurve.obj HighFrequencyAudioCurve.obj ConstantAudioCurve.obj CompoundAudioCurve.obj Profiler.obj Resampler.obj FFT.obj AudioCurveCalculator.obj getopt_long.obj getopt.obj kiss_fftr.obj kiss_fft.obj resample.obj VectorOpsComplex.obj Thread.obj sysutils.obj Allocators.obj
    elseif strcmpi(varargin{1},'vowel_shifter_smb')
        mex('../vowel_shifter.cpp', '-D__LITTLE_ENDIAN__', '-D__WINDOWS_ASIO__', '-D__WINDOWS_MM__', '-c', '-I../stk-4.5.0/include', path_param, '-DSMBSHIFTER');
        mex ../smbPitchShift.cpp -g -c -I../stk-4.5.0/include
        mex ../dywapitchtrack.c -g -c
        mex vowel_shifter.obj smbPitchShift.obj dywapitchtrack.obj Stk.obj RtAudio.obj FM.obj ADSR.obj FileLoop.obj SineWave.obj Wurley.obj Rhodey.obj BeeThree.obj Drummer.obj OnePole.obj FileWvIn.obj FileRead.obj TwoZero.obj asio.obj asiolist.obj asiodrivers.obj iasiothiscallresolver.obj -L../ -lfftw3-3.lib -output vowel_shifter_smb
    elseif strcmpi(varargin{1},'vowel_shifter_cpv')
        mex('../vowel_shifter.cpp', '-D__LITTLE_ENDIAN__', '-D__WINDOWS_ASIO__', '-D__WINDOWS_MM__', '-c', '-I../stk-4.5.0/include', path_param, '-DCPVSHIFTER');
        mex ../cpvPitchShift.cpp -g -c -I../stk-4.5.0/include
        mex ../dywapitchtrack.c -g -c
        mex vowel_shifter.obj cpvPitchShift.obj dywapitchtrack.obj Stk.obj RtAudio.obj FM.obj ADSR.obj FileLoop.obj SineWave.obj Wurley.obj Rhodey.obj BeeThree.obj Drummer.obj OnePole.obj FileWvIn.obj FileRead.obj TwoZero.obj asio.obj asiolist.obj asiodrivers.obj iasiothiscallresolver.obj -L../ -lfftw3-3.lib -lwinmm -lWsock32 -output vowel_shifter_cpv 
    elseif strcmpi(varargin{1},'vowel_shifter_dirac')
        mex ../vowel_shifter_dirac.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -g -c -I../stk-4.5.0/include
        mex ../dywapitchtrack.c -g -c
        mex vowel_shifter_dirac.obj dywapitchtrack.obj Stk.obj RtAudio.obj FM.obj ADSR.obj FileLoop.obj SineWave.obj Wurley.obj Rhodey.obj BeeThree.obj Drummer.obj OnePole.obj FileWvIn.obj FileRead.obj TwoZero.obj asio.obj asiolist.obj asiodrivers.obj iasiothiscallresolver.obj -L../ -lDiracLE.lib
    else
        fprintf('invalid input\n');
    end
    
    cd(currentFolder)
end