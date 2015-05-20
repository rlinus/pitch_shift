#ifndef CPVPITCHSHIFT_H
#define CPVPITCHSHIFT_H

void cpvPitchShiftInit(long stepSize_i, long osamp_i, float sampleRate_i);
void cpvPitchShift(float pitchShift, float *indata, float *outdata);


#endif