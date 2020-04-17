function twoPhase = doublePhase(img_complex)
[sz1, sz2] = size(img_complex);
amp = normMax(abs(img_complex));
cosinv = acos(amp);
phase = angle(img_complex);
twoPhase = zeros(sz1,sz2);
p1 = .5*exp(1j*(phase + cosinv));
p2 = .5*exp(1j*(phase - cosinv));

% create a mosaic from p1 and p2
twoPhase(1:2:end, 1:2:end) = p1(1:2:end, 1:2:end);
twoPhase(2:2:end, 2:2:end) = p1(2:2:end, 2:2:end);
twoPhase(1:2:end, 2:2:end) = p2(1:2:end, 2:2:end);
twoPhase(2:2:end, 1:2:end) = p2(2:2:end, 1:2:end);


