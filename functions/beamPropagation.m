function vol_output = beamPropagation(input, Nstep, propPhase, modulPhase)

[sz1, sz2] = size(propPhase);

if Nstep ~= 0
    vol_output = zeros(sz1, sz2, Nstep);
    if nargin == 4
        vol_output(:,:,1) = input.*modulPhase(:,:,1);
    else
        vol_output(:,:,1) = input;
    end
    for i = 2:Nstep
        propOutput = propImagePhase(input, propPhase);
        input = propOutput;
        if nargin == 4
            vol_output(:,:,i) = propOutput.*modulPhase(:,:,i);
        else
            vol_output(:,:,i) = propOutput;
        end
    end
else
    vol_output = input;
end