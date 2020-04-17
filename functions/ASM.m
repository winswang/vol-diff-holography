function I = ASM(E0, X0, Y0, z, lambda, mask)
% Angular spectrum method
% the function receives a field E0 at wavelength lambda
% and returns the field E after distance z
% dx, dy are spatial resolution
k = 2*pi/lambda;
[row, col] = size(E0);

dx = X0/col; dy = Y0/row; % spatial frequency
Kx = 2*pi/dx; Ky = 2*pi/dy;
kx = linspace((-Kx/2), (Kx/2), col);
ky = linspace((-Ky/2), (Ky/2), row);

[kxgrid, kygrid] = meshgrid(kx, ky);

% construct the circle function
circ = sqrt(kxgrid.^2 + kygrid.^2)/k;
circ(circ>1) = 0;
circ(circ<=1) = 1;

F = fftshift(fft2(ifftshift(E0)));
factor = exp(1i*z*sqrt(k^2 - kxgrid.^2 - kygrid.^2));
E = fftshift(ifft2(ifftshift(F.*factor.*circ)));

if nargin <6
    mask = ones(row, col);
end

I = real(E);
I = I.*mask;
end