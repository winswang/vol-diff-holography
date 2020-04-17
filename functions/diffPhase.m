function diff = diffPhase(phase3D)
[sz1, sz2, sz3] = size(phase3D);
diff = zeros(sz1, sz2, sz3);
for i = 1:(sz3-1)
    diff(:,:,i+1) = phase3D(:,:,i+1) - phase3D(:,:,i);
end