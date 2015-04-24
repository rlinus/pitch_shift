function make(varargin)
    cd(fileparts(which(mfilename)));
    
    if isempty(varargin{1}) || ~ischar(varargin{1})
        varargin{1} = 'all';
    end
    
    if strcmpi(varargin{1},'all')
        % compile extern code
        mex stk-4.5.0/src/include/asio.cpp -DWIN32 -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -Istk-4.5.0/src/include
        mex stk-4.5.0/src/include/asiolist.cpp  -DWIN32 -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -Istk-4.5.0/src/include
        mex stk-4.5.0/src/include/asiodrivers.cpp  -DWIN32 -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -Istk-4.5.0/src/include
        mex stk-4.5.0/src/include/iasiothiscallresolver.cpp  -DWIN32 -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -c -Istk-4.5.0/src/include

        mex stk-4.5.0/src/Stk.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -g -c -Istk-4.5.0/src/include -Istk-4.5.0/include
        mex stk-4.5.0/src/RtAudio.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -g -c -Istk-4.5.0/src/include -Istk-4.5.0/include

        mex rubberband-1.8.1/src/RubberBandStretcher.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1

        mex rubberband-1.8.1/src/StretchCalculator.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1
        mex rubberband-1.8.1/src/StretcherChannelData.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1
        mex rubberband-1.8.1/src/StretcherImpl.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1
        mex rubberband-1.8.1/src/StretcherProcess.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1

        mex rubberband-1.8.1/src/audiocurves/SpectralDifferenceAudioCurve.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src
        mex rubberband-1.8.1/src/audiocurves/SilentAudioCurve.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src
        mex rubberband-1.8.1/src/audiocurves/PercussiveAudioCurve.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src
        mex rubberband-1.8.1/src/audiocurves/HighFrequencyAudioCurve.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src
        mex rubberband-1.8.1/src/audiocurves/ConstantAudioCurve.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src
        mex rubberband-1.8.1/src/audiocurves/CompoundAudioCurve.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src

        mex rubberband-1.8.1/src/base/Profiler.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src

        mex rubberband-1.8.1/src/dsp/Resampler.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src
        mex rubberband-1.8.1/src/dsp/FFT.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src
        mex rubberband-1.8.1/src/dsp/AudioCurveCalculator.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src

        mex rubberband-1.8.1/src/getopt/getopt_long.c -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src
        mex rubberband-1.8.1/src/getopt/getopt.c -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src

        mex rubberband-1.8.1/src/kissfft/kiss_fftr.c -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src
        mex rubberband-1.8.1/src/kissfft/kiss_fft.c -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src

        mex rubberband-1.8.1/src/speex/resample.c -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src

        mex rubberband-1.8.1/src/system/VectorOpsComplex.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src
        mex rubberband-1.8.1/src/system/Thread.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src
        mex rubberband-1.8.1/src/system/sysutils.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src
        mex rubberband-1.8.1/src/system/Allocators.cpp -DWIN32  -D__MSVC__ -DUSE_KISSFFT -DUSE_SPEEX -c -Irubberband-1.8.1 -Irubberband-1.8.1/src
    end
    
    % compile custom code
    mex rt_pitch_shifter.cpp -D__LITTLE_ENDIAN__ -D__WINDOWS_ASIO__ -D__WINDOWS_MM__ -g -c -Istk-4.5.0/include -Irubberband-1.8.1/rubberband

    % link everything
    mex rt_pitch_shifter.obj Stk.obj RtAudio.obj asio.obj asiolist.obj asiodrivers.obj iasiothiscallresolver.obj RubberBandStretcher.obj StretchCalculator.obj StretcherChannelData.obj StretcherImpl.obj StretcherProcess.obj SpectralDifferenceAudioCurve.obj SilentAudioCurve.obj PercussiveAudioCurve.obj HighFrequencyAudioCurve.obj ConstantAudioCurve.obj CompoundAudioCurve.obj Profiler.obj Resampler.obj FFT.obj AudioCurveCalculator.obj getopt_long.obj getopt.obj kiss_fftr.obj kiss_fft.obj resample.obj VectorOpsComplex.obj Thread.obj sysutils.obj Allocators.obj
end