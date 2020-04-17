function ang_phi = rampPhase2D(params, ang_xy, z_dim)
x_axis = linspace(-params.X/2, params.X/2, params.Nxy);
z_axis = z_dim*ones(1, params.Nxy);
[xgrid, ygrid] = meshgrid(x_axis, x_axis);
[xgrid, zgrid] = meshgrid(x_axis, z_axis);
ang_phi = exp(1i*params.k*(xgrid*sin(ang_xy(1))+zgrid*cos(ang_xy(1))));
% ang_phiy = exp(1i*params.k*(ygrid*sin(ang_xy(2))+zgrid*cos(ang_xy(2))));
% ang_phi = ang_phix + ang_phiy;