#ifndef SMBPITCHSHIFT_H
#define SMBPITCHSHIFT_H

void smbPitchShiftInit(long stepSize_i, long osamp_i, double sampleRate_i);
void smbPitchShift(double pitchShift, double *indata, double *outdata);


#endif