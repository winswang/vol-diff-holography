function y3D = shiftZ(x3D)
% input signal is a 3D signal
[sz1, sz2, sz3] = size(x3D);
midIndex = floor(sz3/2);
firstHalf = x3D(:,:,1:midIndex);
secondHalf = x3D(:,:,(midIndex+1):end);
y3D = cat(3, secondHalf, firstHalf);
