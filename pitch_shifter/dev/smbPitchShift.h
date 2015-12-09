/* smbPitchShift.h
 *
 * smbPitchShift is a real-time pitch shifter based on Stephan Bernsee's modified
 * phase vocoder algorithm (http://blogs.zynaptiq.com/bernsee/pitch-shifting-using-the-ft/).
 *
 * Initalize with smbPitchShiftInit() befor using cpvPitchShift(). stepSize_i and fftFrameSize_i
 * must be powers of 2 and 2*stepSize_i <= fftFrameSize_i <= 8192. stepSize_i is the hopsize
 * of the phase vocoder. pitchShift is the pitch shift factor
 * and it must hold: 0.5 <= pitchShift <= 2. indata and outdata must point to arrays of size stepSize_i.
 */


#ifndef SMBPITCHSHIFT_H
#define SMBPITCHSHIFT_H

void smbPitchShiftInit(long stepSize_i, long fftFrameSize_i, double sampleRate_i);
void smbPitchShift(double pitchShift, double *indata, double *outdata);


#endif