%% multiplexing
clc; close all; clear all
addpath('functions');
cm = pink(256);

% define parameters
% length unit in um
params.wavelength = .532; % in um
params.n0 = 1.5; % refractive index of hologram
params.dn_max = 0.05; % max variation for n
params.k = 2*pi/params.wavelength*params.n0; % wave number
params.k_max = 6*params.k/2/pi;
params.Nxy = 256; % number of pixels in x and y
params.dxy = 1/params.k_max;
params.dz = 1/params.k_max;
params.Nz = 256;
params.X = params.Nxy*params.dxy; % total length

% params.dz = .532/5; % depth size in um
% params.Nz = ceil(16/params.dz); % number of depth layers
params.Z = params.dz*params.Nz;
centerZ = floor(params.Nz/2)+1; % index of the central depth layer

% construct Green's Function
G.xyAxis = linspace(-params.X/2, params.X/2, params.Nxy);
G.zAxis = linspace(-params.Z/2, params.Z/2, params.Nz);

[G.xGrid, G.yGrid, G.zGrid] = meshgrid(G.xyAxis, G.xyAxis, G.zAxis);
G.radiusFwd = sqrt(G.xGrid.^2 + G.yGrid.^2 + (G.zGrid + params.Z/2).^2);
G.radiusBwd = sqrt(G.xGrid.^2 + G.yGrid.^2 + (G.zGrid - params.Z/2).^2);

G.fwd = exp(1i*params.k.*G.radiusFwd)./G.radiusFwd;
G.bwd = exp(1i*params.k.*G.radiusBwd)./G.radiusBwd;

G.fwdFFT = fftn(G.fwd);
G.fwdFFTabs = abs(G.fwdFFT);
G.mask_max = 1e4;
G.maskFwdFFT = double(G.fwdFFTabs>G.mask_max);
G.fwdFFT = normSum(G.maskFwdFFT);

G.bwdFFT = fftn(G.bwd);
G.bwdFFTabs = abs(G.bwdFFT);
G.maskBwdFFT = double(G.bwdFFTabs>G.mask_max);
G.bwdFFT = normSum(G.maskBwdFFT);

%% num-8
% recording
ref.xshift = [-4, -3, -2, 2, 3, 4, 5]*5;
% ref.yshift = [-14, 2, 10, 18, 26, -30, -22];
ref.img = ones(params.Nxy, params.Nxy)/(params.Nxy*params.Nxy);
hol.index3D = params.n0*ones(params.Nxy, params.Nxy, params.Nz);
for i = 1:7
    obj.img(:,:,i) = normSum(rgb2gray(im2double(imresize(imread(sprintf('images/num-8-%1d.png',i))...
        , [params.Nxy, params.Nxy]))));
    obj.amp = normSum(sqrt(imgaussfilt(obj.img(:,:,i), 1)));
    obj.input3D = zeros(params.Nxy, params.Nxy, params.Nz);
    obj.input3D(:,:,centerZ) = obj.amp;
    obj.output3DFFT = G.fwdFFT.*fftn(obj.input3D);
    obj.output3D = ifftshift(ifftn(obj.output3DFFT));
    oimg = fftshift(obj.output3D);
    basis.ub(:,:,i) = oimg(:,:,centerZ);
    
    ref.imgFFT = circshift(fft2(ref.img), [ref.xshift(i), 0]);
    ref.ramp = ifft2(ref.imgFFT);
    ref.input3D = zeros(params.Nxy, params.Nxy, params.Nz);
    ref.input3D(:,:,centerZ) = ref.ramp;
    ref.output3DFFT = G.bwdFFT.*fftn(ref.input3D);
    ref.output3D = ifftshift(ifftn(ref.output3DFFT));
    basis.ur(:,:,i) = ref.ramp;
    
    hol.fac = 1;
    hol.intensity3D = complexSquare(obj.output3D + hol.fac*ref.output3D);
    hol.scale = params.dn_max/max(hol.intensity3D(:));
    hol.index3D = hol.index3D + hol.scale*hol.intensity3D;
end

potential3D = 1/4/pi*params.k^2*((hol.index3D/params.n0).^2 - 1);
potential3DFFT = fftshift(fftn(potential3D));
potentialFFTsec = squeeze(potential3DFFT(:,params.Nxy/2+1,:));
figure(101);
imagesc(abs(potentialFFTsec)); axis equal off;
colormap(cm); caxis([0, 2e-3*max(abs(potentialFFTsec(:)))]);

%% playback
code.n1 = [0,0,1,1,0,0,0]/7;
code.n2 = [0,3,1,0,3,3,1]/7;
code.n3 = [0,3,1,1,3,0,1]/7;
code.n4 = [2.8,0,1,1,0,0,1]/7;
code.n5 = [3,3,0,1,3,0,1]/7;
code.n6 = [3,3,0,1,3,3,1]/7;
code.n7 = [0,2.5,1,1,0,0,0]/7;
code.n8 = [3,3,1,1,3,3,1]/7;
code.n9 = [3,3,1,1,3,0,1]/7;
code.n0 = [3,2.5,1,1,2.5,3,0]/7;
code.bar257 = [0,3,0,0,3,0,1]/7;
code.bar1 = [1,0,0,0,0,0,0]/7;
code.bar7 = [0,0,0,0,0,0,1]/7;
pb.code = code.bar257;
pb.codename = 'bar257-adj';
pb.ref = zeros(params.Nxy, params.Nxy);
for i = 1:7
    pb.ref = pb.ref + pb.code(i)*basis.ur(:,:,i);
end

figure(102);
imagesc(angle(pb.ref));

pb.ref3D = zeros(params.Nxy, params.Nxy, params.Nz);
pb.ref3D(:,:,centerZ) = pb.ref;
pb.ref3DFFT = G.bwdFFT.*fftn(pb.ref3D);
pb.ref3D = ifftshift(ifftn(pb.ref3DFFT));
pb.in3D = potential3D.*pb.ref3D;
pb.inFFT = fftshift(fftn(pb.in3D));
pb.scatter3D = ifftshift(ifftn(G.fwdFFT.*fftn(pb.in3D)));
pb.s = squeeze(pb.scatter3D(:,:,centerZ));
figure(103)
imagesc(abs(squeeze(pb.inFFT(:,params.Nxy/2+1,:))));
caxis([0,2e-3*max(abs(pb.inFFT(:)))]); colormap(cm)
figure(104);
imagesc(complexSquare(pb.s)); colormap(cm)

%% save image
imwrite(normMax(angle(pb.ref))*255, cm, sprintf('results/multi-diff-holo-3-30/pb-%s-ref-ang.png', pb.codename));
imwrite(normMax(abs(squeeze(pb.inFFT(:,params.Nxy/2+1,:))))*5e2*255, cm, sprintf('results/multi-diff-holo-3-30/pb-%s-vr-abs.png', pb.codename));
imwrite(normMax(complexSquare(pb.s))*255, cm, sprintf('results/multi-diff-holo-3-30/pb-%s-csq.png', pb.codename));

%%


holo.interfere3D = ifftshift(ifftn(holo.interfere3DFFT));
holo.intensity3D = complexSquare(holo.interfere3D);
holo.scale = params.dn_max/max(holo.intensity3D(:));
holo.index3D = params.n0 + holo.scale*holo.intensity3D;
% separate terms
index3DDC = params.n0 + holo.scale*(complexSquare(obj.output3D) + complexSquare(fac*ref.output3D));
index3DORc = params.n0 + holo.scale*(obj.output3D.*conj(fac*ref.output3D));
index3DOcR = params.n0 + holo.scale*(conj(obj.output3D).*fac.*ref.output3D);

holo.potential3D = 1/4/pi*params.k^2*((holo.index3D/params.n0).^2 - 1);
% separate terms
potentialDC = 1/4/pi*params.k^2*((index3DDC/params.n0).^2 - 1);
potentialORc = 1/4/pi*params.k^2*((index3DORc/params.n0).^2 - 1);
potentialOcR = 1/4/pi*params.k^2*((index3DOcR/params.n0).^2 - 1);

VFFT = fftshift(fftn(holo.potential3D));
v_max = max(abs(VFFT(:)));
VdcFFT = fftshift(fftn(potentialDC));
VorcFFT = fftshift(fftn(potentialORc));
VocrFFT = fftshift(fftn(potentialOcR));
VsumFFT = VdcFFT + VorcFFT + VocrFFT;
figure(103); % potential FFT
subplot(221);
imagesc(abs(squeeze(VFFT(:,params.Nxy/2+1,:)))); axis equal off;
colormap(cm); title('|F\{V\}|'); caxis([0, 1e-2*v_max]);
subplot(222);
imagesc(abs(squeeze(VdcFFT(:,params.Nxy/2+1,:)))); axis equal off;
colormap(cm); title('|F\{V_{DC}\}|'); caxis([0, 1e-2*v_max]);
subplot(223);
imagesc(abs(squeeze(VsumFFT(:,params.Nxy/2+1,:)))); axis equal off;
colormap(cm); title('|F\{V_{sum}\}|'); caxis([0, 1e-2*v_max]);
subplot(224);
imagesc(abs(squeeze(VorcFFT(:,params.Nxy/2+1,:)))); axis equal off;
colormap(cm); title('|F\{V_{O.R*}\}|'); caxis([0, 1e-2*v_max]);

% playback
pb.incident3D = holo.potential3D.*ref.output3D;
pb.incident3DFFT = fftshift(fftn(pb.incident3D));
incidentDCFFT = fftshift(fftn(potentialDC.*ref.output3D));
incidentORcFFT = fftshift(fftn(potentialORc.*ref.output3D));
incidentOcRFFT = fftshift(fftn(potentialOcR.*ref.output3D));

figure(104);
c_max104 = max(max(abs(squeeze(pb.incident3DFFT(:,params.Nxy/2+1,:)))));
subplot(221);
imagesc(abs(squeeze(pb.incident3DFFT(:,params.Nxy/2+1,:)))); axis equal off
title('|F\{V.R\}|'); caxis([0,1e-2*c_max104]);
subplot(222);
imagesc(abs(squeeze(incidentDCFFT(:,params.Nxy/2+1,:)))); axis equal off
title('|F\{V_{DC}.R\}|'); caxis([0,1e-2*c_max104]);
subplot(223);
imagesc(abs(squeeze(incidentORcFFT(:,params.Nxy/2+1,:)))); axis equal off
title('|F\{V_{O.Rc}.R\}|'); caxis([0,1e-2*c_max104]);
subplot(224);
imagesc(abs(squeeze(incidentOcRFFT(:,params.Nxy/2+1,:)))); axis equal off
title('|F\{V_{Oc.R}.R\}|'); caxis([0,1e-2*c_max104]);

pb.scatter3D = ifftshift(ifftn(G.fwdFFT.*fftshift(pb.incident3DFFT)));
scatterDC = ifftshift(ifftn(G.fwdFFT.*fftshift(incidentDCFFT)));
scatterORc = ifftshift(ifftn(G.fwdFFT.*fftshift(incidentORcFFT)));
scatterOcR = ifftshift(ifftn(G.fwdFFT.*fftshift(incidentOcRFFT)));
% pb.total3D = ref.output3D + pb.scatter3D;

figure(115);
imagesc(abs(squeeze(scatterORc(:,params.Nxy/2+1,:)))); axis equal off
figure(105);
subplot(221);
imagesc(complexSquare(pb.scatter3D(:,:,centerZ))); axis equal off
title('scatter')
subplot(222);
imagesc(abs(scatterDC(:,:,centerZ))); axis equal off
title('DC')
subplot(223);
imagesc(abs(scatterORc(:,:,centerZ))); axis equal off
title('O.Rc')
subplot(224);
imagesc(abs(scatterOcR(:,:,centerZ))); axis equal off
title('Oc.R')
% figure(301)
% imagesc(complexSquare(pb.scatter3D(:,:,centerZ+1)));
%% save images
% imwrite(GFFTsec*255, cm, sprintf('results/diff-holo-3-29/green-fft-fwd.png'));
% imwrite(normMax(complexSquare(o_img))*255, cm, sprintf('results/diff-holo-3-29/obj-csq.png'));
% imwrite(normMax(abs(o_img))*255, cm, sprintf('results/diff-holo-3-29/obj-abs.png'));
% imwrite(normMax(complexSquare(obj.amp))*255, cm, sprintf('results/diff-holo-3-29/obj-original-csq.png'));
% imwrite(normMax(complexSquare(obj.section))*255, cm, sprintf('results/diff-holo-3-29/obj-sec-csq.png'));
% imwrite(normMax(abs(H3DFFTsec))*255, cm, sprintf('results/diff-holo-3-29/O-R-fft-sec-abs.png'));
% imwrite(normMax(abs(squeeze(VFFT(:,params.Nxy/2+1,:))))*1e2*255, cm, sprintf('results/diff-holo-3-29/v-abs.png'));
% imwrite(normMax(complexSquare(scatterOcR(:,:,centerZ)))*255, cm, sprintf('results/diff-holo-3-29/s-oc-r-csq.png'));

