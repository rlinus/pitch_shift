#ifndef SMBPITCHSHIFT_H
#define SMBPITCHSHIFT_H

void smbFft(float *fftBuffer, long fftFrameSize, long sign);
double smbAtan2(double x, double y);
void smbPitchShift(float pitchShift, long numSampsToProcess, long fftFrameSize, long osamp, float sampleRate, float *indata, float *outdata);
void smbPitchShift2(float pitchShift, long stepSize, long osamp, float sampleRate, float *indata, float *outdata);


#endif