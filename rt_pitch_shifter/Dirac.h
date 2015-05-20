/*--------------------------------------------------------------------------------
 
 Dirac.h
 Version 3.6.2 [1305130900]
 
 Copyright (C) 2005-2013 The DSP Dimension,
 Stephan M. Bernsee (SMB)
 All rights reserved
 
 CONFIDENTIAL: This document contains confidential information. 
 Do not disclose any information contained in this document to any
 third-party without the prior written consent of Stephan M. Bernsee/
 The DSP Dimension.
  
 --------------------------------------------------------------------------------*/

// This file contains all the constants and prototypes for the "real" DIRAC calls
// that you will need in your project.


#ifndef __DIRAC__
#define __DIRAC__



// Windows DLL definitions
// ----------------------------------------------------------------------------

#ifndef __APPLE__
	#ifdef DIRAC_AS_DLL
		#define DLL_DEF_TYPE __declspec(dllexport)
	#else
		#define DLL_DEF_TYPE
	#endif
#else
	#define DLL_DEF_TYPE __attribute__((visibility("default")))
#endif


// Function prototypes
// ----------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif

	/* ******************* DIRAC CORE API ********************** */
	/* Dirac Core calls */
	DLL_DEF_TYPE void *DiracCreate(long lambda, long quality, long numChannels, float sampleRateHz, long (*readFromChannelsCallback)(float **data, long numFrames, void *userData), void *userData);
	DLL_DEF_TYPE void *DiracCreateInterleaved(long lambda, long quality, long numChannels, float sampleRateHz, long (*readFromInterleavedChannelsCallback)(float *data, long numFrames, void *userData), void *userData);
	DLL_DEF_TYPE long DiracSetProperty(long selector, long double value, void *dirac);
	DLL_DEF_TYPE long double DiracGetProperty(long selector, void *dirac);
	DLL_DEF_TYPE void DiracReset(bool clear, void *dirac);
	DLL_DEF_TYPE long DiracProcess(float **audioOut, long numFrames, void *dirac);
	DLL_DEF_TYPE long DiracProcessInterleaved(float *audioOut, long numFrames, void *dirac);
	DLL_DEF_TYPE void DiracDestroy(void *dirac);
	DLL_DEF_TYPE void DiracSetProcessingBeganCallback(void (*processingCallback)(unsigned long position, void *userData), void *userData, void *dirac);
	
	// available in Dirac PRO only	
	DLL_DEF_TYPE long DiracSetTuningTable(float *frequencyTable, long numFrequencies, void *dirac);

	
	
	/* ******************* DIRAC RETUNE API ********************** */
	/* This is the interface for DIRAC3's Retune algorithm */
	DLL_DEF_TYPE void *DiracRetuneCreate(long quality, float sampleRateHz, float referenceTuningHz);
	DLL_DEF_TYPE void DiracRetuneDestroy(void *instance);
	DLL_DEF_TYPE void DiracRetuneProcess(short *indata, short *outdata, long numSampsToProcess, void *instance);
	DLL_DEF_TYPE void DiracRetuneProcessFloat(float *indata, float *outdata, long numSampsToProcess, void *instance);
	DLL_DEF_TYPE void DiracRetuneSetKeyList(float *tuningCentRelativeToKey0, long numKeysPerOctave, long octaveOffsetKeyNo, void *diracRetune);
	DLL_DEF_TYPE void DiracRetuneSetProperties(float correctionAmountPercent,
								  float correctionCaptureCent,
								  float correctionAutoBypassThreshold,
								  float correctionAmbienceThreshold,
								  void *instance);
	
	DLL_DEF_TYPE float DiracRetuneGetPitchHz(void *instance);
	DLL_DEF_TYPE bool DiracRetuneGetKeyStatus(long keyNo, void *instance);
	DLL_DEF_TYPE void DiracRetuneSetKeyStatus(long keyNo, bool enable, void *instance);
	DLL_DEF_TYPE unsigned long DiracRetuneGetAllowedKeysMask(void *instance);
	DLL_DEF_TYPE void DiracRetuneSetAllowedKeysMask(unsigned long mask, void *instance);
	DLL_DEF_TYPE float DiracRetuneGetClosestKeyDetuneCent(bool respectKeyState, void *instance);
	DLL_DEF_TYPE long DiracRetuneGetClosestKey(bool respectKeyState, void *instance);
	DLL_DEF_TYPE void DiracRetunePrintInternalTuningTable(void *instance);
	DLL_DEF_TYPE long DiracRetuneLatencyFrames(float sampleRateHz);

	/* Deprecated calls */
	DLL_DEF_TYPE void DiracRetuneSetPitchHz(float pitchHz, void *instance);
	DLL_DEF_TYPE void DiracRetuneSetTuningReferenceHz(float referenceTuningHz, void *instance);
	DLL_DEF_TYPE void DiracRetuneSetTuningTable(float *frequencyTable, long numFrequencies, void *instance);
	
	
	
	
	/* ******************* DIRAC FX API ********************** */
	/* This is the interface for DIRAC 3.5's FX mode */
	DLL_DEF_TYPE void *DiracFxCreate(long quality, float sampleRateHz, long numChannels);
	DLL_DEF_TYPE long DiracFxMaxOutputBufferFramesRequired(long double timeFactor, long double pitchFactor, long numInputFrames);
	DLL_DEF_TYPE long DiracFxOutputBufferFramesRequiredNextCall(long double timeFactor, long double pitchFactor, long numInputFrames, void *instance);
	DLL_DEF_TYPE long DiracFxLatencyFrames(float sampleRateHz);
	DLL_DEF_TYPE void DiracFxDestroy(void *instance);
	DLL_DEF_TYPE long DiracFxProcessFloat(long double timeFactor, long double pitchFactor, float **indata, float **outdata, long numInputFrames, void *instance);
	DLL_DEF_TYPE long DiracFxProcessFloatInterleaved(long double timeFactor, long double pitchFactor, float *indata, float *outdata, long numInputFrames, void *instance);
	DLL_DEF_TYPE long DiracFxProcess(long double timeFactor, long double pitchFactor, short **indata, short **outdata, long numInputFrames, void *instance);
	DLL_DEF_TYPE long DiracFxProcessInterleaved(long double timeFactor, long double pitchFactor, short *indata, short *outdata, long numInputFrames, void *instance);
	DLL_DEF_TYPE void DiracFxReset(bool clear, void *instance);
	

	
	/* Utilities */
	DLL_DEF_TYPE const char *DiracVersion(void);
	DLL_DEF_TYPE void DiracStartClock(void);
	DLL_DEF_TYPE long double DiracClockTimeSeconds(void);
	DLL_DEF_TYPE float DiracPeakCpuUsagePercent(void *dirac);
	DLL_DEF_TYPE long double DiracValidateStretchFactor(long double factor);
	DLL_DEF_TYPE void DiracPrintSettings(void *dirac);
	DLL_DEF_TYPE const char *DiracErrorToString(long error);

	
#ifdef __cplusplus
}
#endif



// Property enums
// ----------------------------------------------------------------------------

enum
{
	kDiracPropertyPitchFactor = 100,
	kDiracPropertyTimeFactor,
	kDiracPropertyFormantFactor,
	kDiracPropertyCompactSupport,
	kDiracPropertyCacheGranularity,
	kDiracPropertyCacheMaxSizeFrames,
	kDiracPropertyCacheNumFramesLeftInCache,
	kDiracPropertyUseConstantCpuPitchShift,
	kDiracPropertyDoPitchCorrection,
	kDiracPropertyOutputGainDb,
	kDiracPropertyPitchCorrectionBasicTuningHz = 400,
	kDiracPropertyPitchCorrectionSlurTime,
	kDiracPropertyPitchCorrectionDoFormantCorrection = 500,
	kDiracPropertyPitchCorrectionFundamentalFrequency,
	
	kDiracPropertyNumProperties
};


// Lambda enums
// ----------------------------------------------------------------------------
enum
{
	kDiracLambdaPreview = 200,
	kDiracLambda1,
	kDiracLambda2,
	kDiracLambda3,
	kDiracLambda4,
	kDiracLambda5,
	kDiracLambdaTranscribe,
	
	kDiracPropertyNumLambdas
};




// Quality enums
// ----------------------------------------------------------------------------
enum
{
	kDiracQualityPreview = 300,	
	kDiracQualityGood,
	kDiracQualityBetter,
	kDiracQualityBest,
	
	kDiracPropertyNumQualities
};




// Error enums
// ----------------------------------------------------------------------------
enum
{
	kDiracErrorNoErr		= 0,	
	kDiracErrorParamErr		= -1,
	kDiracErrorUnknownErr	= -2,
	kDiracErrorInvalidCb	= -3,
	kDiracErrorCacheErr		= -4,
	kDiracErrorNotInited	= -5,
	kDiracErrorMultipleInits	= -6,
	kDiracErrorFeatureNotSupported	= -7,
	kDiracErrorMemErr		= -108,
	kDiracErrorDemoTimeoutReached = -10001,
	
	kDiracErrorNumErrs
};





#endif /* __DIRAC__ */


