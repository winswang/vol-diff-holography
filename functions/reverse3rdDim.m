function new_x = reverse3rdDim(x)

[sz1, sz2, sz3] = size(x);

new_x = zeros(sz1, sz2, sz3);
for i = 1:sz3
    new_x(:,:,i) = x(:,:,sz3-i+1);
end