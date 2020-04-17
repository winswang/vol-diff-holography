function phi = constructPhaseKernel(params, dz)

k = params.k;
Kxy = 2*pi/params.dxy;
kxy = linspace((-Kxy/2), (Kxy/2), params.Nxy);
[kxgrid, kygrid] = meshgrid(kxy, kxy);
% construct the circle function
circ = sqrt(kxgrid.^2 + kygrid.^2)/k;
circ(circ>1) = 0;
circ(circ<=1) = 1;

phi = exp(1i*dz*sqrt(k^2 - kxgrid.^2 - kygrid.^2)).*circ;