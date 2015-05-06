#ifndef SMBPITCHSHIFT_H
#define SMBPITCHSHIFT_H

void smbFft(float *fftBuffer, long fftFrameSize, long sign);
double smbAtan2(double x, double y);
void smbPitchShiftInit(long stepSize_i, long osamp_i, float sampleRate_i);
void smbPitchShift(float pitchShift, float *indata, float *outdata);


#endif