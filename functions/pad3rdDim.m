function xout = pad3rdDim(x, val)
[s1, s2, s3] = size(x);
if nargin == 1
    val = 0;
end

padCube = ones(s1, s2, s3)*val;

xout = cat(3, padCube, x, padCube);
