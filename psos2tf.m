function [B,A] = psos2tf(sos,g)

if nargin<2, g=1; end

[nsecs,tmp] = size(sos);
if nsecs<1, B=[]; A=[]; return; end
Bs = sos(:,1:3);
As = sos(:,4:6);
B = Bs(1,:);
A = As(1,:);
for i=2:nsecs
  B = conv(B,As(i,:)) + conv(A,Bs(i,:));
  A = conv(A,As(i,:));
end