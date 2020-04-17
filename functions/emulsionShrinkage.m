function index3D = emulsionShrinkage(index3D, shrink)
if nargin < 2
    shrink = 0.1; % percentage
end
[sz1, sz2, sz3] = size(index3D);
shk_oneside = floor(shrink/2*sz3);
shk_sz = [sz1,sz3 - 2*shk_oneside];
for i = 1:sz2
    slice = squeeze(index3D(:,i,:));
    shk = imresize(slice, shk_sz);
    shkpad = padarray(shk, [0, shk_oneside], 1, 'both');
    index3D(:,i,:) = shkpad;
end