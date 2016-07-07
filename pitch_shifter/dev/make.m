function make(varargin)
    %this script can be used to compile and link the pitch shifter
    %functions


    currentFolder = pwd;
    cd(fileparts(which(mfilename)));
    mfilefolder = pwd;
    
    if ~exist('bin','dir'), mkdir('bin'); end;
    cd([mfilefolder '/bin']);
    
    if isempty(varargin{1}) || ~ischar(varargin{1})
        fprintf('specify target');
        return;
    end
    
    stk_ver = '4.5.1';
    stk_path = ['../stk-' stk_ver '/src/'];
    stk_include_path = ['-I../stk-' stk_ver '/include'];
    stk_asio_include_path = ['-I../stk-' stk_ver '/src/include'];
    
    if strcmpi(varargin{1},'libs')
        % compile extern code
        
        stk_asio_files = dir([stk_path '/include/*.cpp']);
        for i=1:length(stk_asio_files)
            mex([stk_path '/include/' stk_asio_files(i).name],'-DWIN32', '-D__LITTLE_ENDIAN__', '-D__WINDOWS_ASIO__', '-D__WINDOWS_MM__', '-c', stk_asio_include_path);
        end
        
        stk_files = dir([stk_path '/*.cpp']);
        for i=1:length(stk_files)
            mex([stk_path '/' stk_files(i).name],'-D__LITTLE_ENDIAN__', '-D__WINDOWS_ASIO__', '-D__WINDOWS_MM__', '-c', stk_include_path, stk_asio_include_path);
        end

        %mex([stk_path '/include/*.cpp'],'-DWIN32', '-D__LITTLE_ENDIAN__', '-D__WINDOWS_ASIO__', '-D__WINDOWS_MM__', '-c', stk_asio_include_path);
        %mex([stk_path '/*.cpp'],'-D__LITTLE_ENDIAN__', '-D__WINDOWS_ASIO__', '-D__WINDOWS_MM__', '-c', stk_include_path, stk_asio_include_path);
        
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
    elseif strcmpi(varargin{1},'PsychPitchShifter')
        mex('../PsychPitchShifter.cpp', '-D__LITTLE_ENDIAN__', '-D__WINDOWS_ASIO__', '-D__WINDOWS_MM__', '-c', stk_include_path, '-I../rubberband-1.8.1/rubberband');
        mex ../smbPitchShift.cpp -c 
        mex ../dywapitchtrack.c -c
        mex ../cpvPitchShift.cpp -c
        mex PsychPitchShifter.obj smbPitchShift.obj cpvPitchShift.obj dywapitchtrack.obj Stk.obj RtAudio.obj FM.obj ADSR.obj FileLoop.obj Delay.obj DelayL.obj SineWave.obj Wurley.obj Rhodey.obj BeeThree.obj Drummer.obj OnePole.obj LentPitShift.obj PitShift.obj FileWvIn.obj FileRead.obj TwoZero.obj asio.obj asiolist.obj asiodrivers.obj iasiothiscallresolver.obj RubberBandStretcher.obj StretchCalculator.obj StretcherChannelData.obj StretcherImpl.obj StretcherProcess.obj SpectralDifferenceAudioCurve.obj SilentAudioCurve.obj PercussiveAudioCurve.obj HighFrequencyAudioCurve.obj ConstantAudioCurve.obj CompoundAudioCurve.obj Profiler.obj Resampler.obj FFT.obj AudioCurveCalculator.obj getopt_long.obj getopt.obj kiss_fftr.obj kiss_fft.obj resample.obj VectorOpsComplex.obj Thread.obj sysutils.obj Allocators.obj -L../ -lfftw3-3.lib -output ../../PsychPitchShifter

        
        %mex('../PsychTest.cpp', '-D__LITTLE_ENDIAN__', '-D__WINDOWS_ASIO__', '-D__WINDOWS_MM__', '-c', stk_include_path, '-I../rubberband-1.8.1/rubberband');
        %mex PsychTest.obj smbPitchShift.obj cpvPitchShift.obj dywapitchtrack.obj Stk.obj RtAudio.obj FM.obj ADSR.obj FileLoop.obj Delay.obj DelayL.obj SineWave.obj Wurley.obj Rhodey.obj BeeThree.obj Drummer.obj OnePole.obj LentPitShift.obj PitShift.obj FileWvIn.obj FileRead.obj TwoZero.obj asio.obj asiolist.obj asiodrivers.obj iasiothiscallresolver.obj RubberBandStretcher.obj StretchCalculator.obj StretcherChannelData.obj StretcherImpl.obj StretcherProcess.obj SpectralDifferenceAudioCurve.obj SilentAudioCurve.obj PercussiveAudioCurve.obj HighFrequencyAudioCurve.obj ConstantAudioCurve.obj CompoundAudioCurve.obj Profiler.obj Resampler.obj FFT.obj AudioCurveCalculator.obj getopt_long.obj getopt.obj kiss_fftr.obj kiss_fft.obj resample.obj VectorOpsComplex.obj Thread.obj sysutils.obj Allocators.obj -L../ -lfftw3-3.lib -output ../../PsychTest
    elseif strcmpi(varargin{1},'PitchLoop')
        mex('../PitchLoop.cpp', '-D__LITTLE_ENDIAN__', '-D__WINDOWS_ASIO__', '-D__WINDOWS_MM__', '-c', stk_include_path, '-I../rubberband-1.8.1/rubberband');
        mex ../smbPitchShift.cpp -c
        mex ../cpvPitchShift.cpp -c
        mex PitchLoop.obj smbPitchShift.obj cpvPitchShift.obj Stk.obj RtAudio.obj FM.obj ADSR.obj FileLoop.obj Delay.obj DelayL.obj SineWave.obj Wurley.obj Rhodey.obj BeeThree.obj Drummer.obj OnePole.obj LentPitShift.obj PitShift.obj FileWvIn.obj FileRead.obj TwoZero.obj asio.obj asiolist.obj asiodrivers.obj iasiothiscallresolver.obj RubberBandStretcher.obj StretchCalculator.obj StretcherChannelData.obj StretcherImpl.obj StretcherProcess.obj SpectralDifferenceAudioCurve.obj SilentAudioCurve.obj PercussiveAudioCurve.obj HighFrequencyAudioCurve.obj ConstantAudioCurve.obj CompoundAudioCurve.obj Profiler.obj Resampler.obj FFT.obj AudioCurveCalculator.obj getopt_long.obj getopt.obj kiss_fftr.obj kiss_fft.obj resample.obj VectorOpsComplex.obj Thread.obj sysutils.obj Allocators.obj -L../ -lfftw3-3.lib -output ../../PitchLoop
    elseif strcmpi(varargin{1},'PrintDeviceList')
        mex('../PrintDeviceList.cpp', '-c', stk_include_path);
        mex PrintDeviceList.obj Stk.obj RtAudio.obj asio.obj asiolist.obj asiodrivers.obj iasiothiscallresolver.obj -output ../../PrintDeviceList
    else
        fprintf('invalid input\n');
    end
    
    cd(currentFolder)
end