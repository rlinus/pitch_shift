/* cpvPitchShift.h
 *
 * cpvPitchShift is a real-time pitch shifter based on a traditional phase vocoder.
 *
 * Initalize with cpvPitchShiftInit() befor using cpvPitchShift(). stepSize_i and fftFrameSize_i
 * must be powers of 2 and 2*stepSize_i <= fftFrameSize_i <= 8192. stepSize_i is the hopsize
 * of the phase vocoder. pitchShift is the pitch shift factor
 * and it must hold: 0.5 <= pitchShift <= 2. indata and outdata must point to arrays of size stepSize_i.
 */

#ifndef CPVPITCHSHIFT_H
#define CPVPITCHSHIFT_H

void cpvPitchShiftInit(long stepSize_i, long fftFrameSize_i, double sampleRate_i);
void cpvPitchShift(double pitchShift, double *indata, double *outdata);


#endif