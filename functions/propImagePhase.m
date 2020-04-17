function obj_d = propImagePhase(img, phi)
obj_fft = fftshift(fft2(ifftshift(img)));
obj_d = fftshift(ifft2(ifftshift(obj_fft.*phi)));