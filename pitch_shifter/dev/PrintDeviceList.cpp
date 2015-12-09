#include <iostream>
#include <string>
#include "mex.h"
#include "RtAudio.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    if(nlhs>0) {
            mexErrMsgIdAndTxt("rt_pitch_shifter:nlhs", "Too many output arguments.");
            return;
    }
    
    RtAudio d;
    
    int n = d.getDeviceCount();
    
    if(n==0) mexPrintf("No ASIO devices found!\n");
    
    for(int i=0; i<n;++i){
        RtAudio::DeviceInfo dI = d.getDeviceInfo(i);
        
        mexPrintf("Device %i:\n\tName: %s\n\tProbed: %i\n\tNumber of input channels: %i\n\tNumber of output channels: %i\n\tNumber of duplex channels: %i\n",i,dI.name,dI.probed,dI.outputChannels,dI.inputChannels,dI.duplexChannels);
    }
    return;
}