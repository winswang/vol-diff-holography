function [field3DPad, half_sz1, half_sz2, half_sz3] = pad3D(field3D, mode)
[sz1, sz2, sz3] = size(field3D);
half_sz3 = ceil(sz3/2);
half_sz1 = ceil(sz1/2);
half_sz2 = ceil(sz2/2);
if nargin < 2
    mode = 3;
end
if mode == 3 % only 3rd dimension
    field3DPad = padarray(field3D, [0 0 half_sz3], 0, 'both');
else
    field3DPad = padarray(field3D, [half_sz1 half_sz2 half_sz3], 0, 'both');
end