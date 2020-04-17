function hogel = phase2index(hogel, params)
% a differential operation on depth is applied
hogel.index3D = zeros(params.Nxy, params.Nxy, hogel.Nz);
for i = 1:(hogel.Nz-1)
    hogel.index3D(:,:,i+1) = (hogel.phase3D(:,:,i+1) - hogel.phase3D(:,:,i))...
        /params.k/hogel.dz;
end

hogel.index3D = hogel.n0 + hogel.index3D;