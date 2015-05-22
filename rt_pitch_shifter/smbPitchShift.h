#ifndef SMBPITCHSHIFT_H
#define SMBPITCHSHIFT_H

void smbFft(double *fftBuffer, long fftFrameSize, long sign);
double smbAtan2(double x, double y);
void smbPitchShiftInit(long stepSize_i, long osamp_i, double sampleRate_i);
void smbPitchShift(double pitchShift, double *indata, double *outdata);


#endif