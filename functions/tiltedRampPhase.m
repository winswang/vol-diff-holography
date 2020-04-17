function ang_phi = tiltedRampPhase(params, ang, z_dim)
x_axis = linspace(-params.X/2, params.X/2, params.Nxy);
z_axis = z_dim*ones(1, params.Nxy);
[xgrid, zgrid] = meshgrid(x_axis, z_axis);
ang_phi = exp(1i*params.k*(xgrid*sin(ang)+zgrid*cos(ang)));