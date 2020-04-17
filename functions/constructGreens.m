function params = constructGreens(params)
% params should have:
% params.Nxy
% params.dxy
% params.Nz
% params.dz
% params.k

params.kxy = linspace(-params.dxy*params.Nxy/2, params.dxy*params.Nxy/2, params.Nxy);
% params.kxy = linspace(0, params.dxy*params.Nxy, params.Nxy);
params.kz = linspace(-params.dz*params.Nz/2, params.dz*params.Nz/2, params.Nz);
% params.kzB = linspace(-params.dz*params.Nz, 0, params.Nz);
% params.kzF = linspace(0, params.dz*params.Nz, params.Nz);

[Kx, Ky, Kz] = meshgrid(params.kxy, params.kxy, params.kz);
KzF = Kz + params.dz*params.Nz/2;
KzB = KzF - params.dz*params.Nz;
RF = sqrt(Kx.^2 + Ky.^2 + KzF.^2);
RB = sqrt(Kx.^2 + Ky.^2 + KzB.^2);
R = sqrt(Kx.^2 + Ky.^2 + Kz.^2);

params.fGreens = exp(1i*params.k.*RF)./RF;
params.bGreens = exp(1i*params.k.*RB)./RB;
params.Greens = exp(1i*params.k.*R)./R;
if params.padMode == 3
    params.zeroZ = zeros(params.Nxy, params.Nxy, params.padNz);
    greens = cat(3, params.zeroZ, params.greens, params.zeroZ);
elseif params.padMode == 123
    greens = zeros(params.Nxy+2*params.padNxy, params.Nxy+2*params.padNxy, ...
        params.Nz+2*params.padNz);
    for i = (params.padNz+1):(params.padNz+1+params.Nz)
        greens((params.padNxy+1):(params.padNxy+1+params.Nxy),...
            (params.padNxy+1):(params.padNxy+1+params.Nxy), i) = ...
            params.greens(:,:,i-params.padNz);
    end
else
    fGreens = params.fGreens;
    bGreens = params.bGreens;
end

params.fGreensFFT = fftn(fGreens);
params.bGreensFFT = fftn(bGreens);
params.GreensFFT = fftn(params.Greens);

