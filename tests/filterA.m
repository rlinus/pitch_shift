%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  A-weighting Filter                  %
%              with MATLAB Implementation              %
%                                                      %
% Author: M.Sc. Eng. Hristo Zhivomirov        06/01/14 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xA = filterA(x, fs)

% function: xA = filterA(x, fs)
% x - signal in the time domain
% fs - sampling frequency, Hz
% xA - filtered signal in the time domain

% determine the signal size
sz = size(x);

% represent x as column-vector
x = x(:);

% signal length
xlen = length(x);

% number of unique points
NumUniquePts = ceil((xlen+1)/2);

% FFT
X = fft(x);

% fft is symmetric, throw away second half
X = X(1:NumUniquePts);

% frequency vector with NumUniquePts points
f = (0:NumUniquePts-1)*fs/xlen;

% A-weighting filter coefficients
c1 = 3.5041384e16;
c2 = 20.598997^2;
c3 = 107.65265^2;
c4 = 737.86223^2;
c5 = 12194.217^2;

% evaluate the A-weighting filter in the frequency domain
f = f.^2;
num = c1*f.^4;
den = ((c2+f).^2) .* (c3+f) .* (c4+f) .* ((c5+f).^2);
A = num./den;
A = A(:);

% filtering in the frequency domain
XA = X.*A;

% reconstruct the whole spectrum
if rem(xlen, 2)                     % odd xlen excludes the Nyquist point
    XA = [XA; conj(XA(end:-1:2))];
else                                % even xlen includes the Nyquist point
    XA = [XA; conj(XA(end-1:-1:2))];
end

% IFFT
xA = real(ifft(XA));
    
% represent the filtered signal in the form of the original one
xA = reshape(xA, sz);

end