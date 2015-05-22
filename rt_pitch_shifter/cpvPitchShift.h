#ifndef CPVPITCHSHIFT_H
#define CPVPITCHSHIFT_H

void cpvPitchShiftInit(long stepSize_i, long osamp_i, double sampleRate_i);
void cpvPitchShift(double pitchShift, double *indata, double *outdata);


#endif