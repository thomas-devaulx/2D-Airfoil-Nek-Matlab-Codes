%% Nek5000 Post-Processing: Combined Vorticity + Mesh Overlay
clear; clc; close all;

% Load Data
[data] = readnek('N0012_R5e3_A080.f00314');

X = data(:,:,1); xg = X(:);
Y = data(:,:,2); yg = Y(:);
Z = data(:,:,3); zg = Z(:);
u = data(:,:,4); ug = u(:);
v = data(:,:,5); vg = v(:);
w = data(:,:,6); wg = w(:);
p = data(:,:,7); pg = p(:);

[n_elem, pts_per_elem] = size(X);
N = sqrt(pts_per_elem); 

% Create Interpolants
F_u = scatteredInterpolant(xg, yg, zg, ug, 'linear', 'none');
F_v = scatteredInterpolant(xg, yg, zg, vg, 'linear', 'none');
F_w = scatteredInterpolant(xg, yg, zg, wg, 'linear', 'none');
F_p = scatteredInterpolant(xg, yg, zg, pg, 'linear', 'none');

% Visualization Parameters
clim_val = 2; % limit of the vorticity field colorbar
v_levels = linspace(-clim_val, clim_val, 120); % number of levels, (...,...,60)
cmap = b2r(-clim_val, clim_val); 

mesh_color = [0 0 0 0.05]; 
mesh_width = 0.25;

%% calculate vorticity
% Define 3D evaluation domain
nx = 400; ny = 200; nz = 100;
xi = linspace(-0.5, 5.0, nx);
yi = linspace(-0.5, 0.5, ny);
zi = linspace(min(zg), max(zg), nz);
[XI, YI, ZI] = meshgrid(xi, yi, zi);

% Interpolate velocity components onto 3D grid
UQ = F_u(XI, YI, ZI);
VQ = F_v(XI, YI, ZI);
WQ = F_w(XI, YI, ZI);

% Compute 3D spatial velocity gradients
[dUdx, dUdy, dUdz] = gradient(UQ, xi, yi, zi);
[dVdx, dVdy, dVdz] = gradient(VQ, xi, yi, zi);
[dWdx, dWdy, dWdz] = gradient(WQ, xi, yi, zi);

% Vorticity vector components
vort_x = dWdy - dVdz;
vort_y = dUdz - dWdx;
vort_z = dVdx - dUdy;

% Vorticity magnitude
vort_mag = sqrt(vort_x.^2 + vort_y.^2 + vort_z.^2);

%% Q-criterion
% Calculate Q-criterion 
Q = 0.5 * (vort_x.^2 + vort_y.^2 + vort_z.^2) - ...
    (dUdx.^2 + dVdy.^2 + dWdz.^2 + 2*(dUdy.*dVdx + dUdz.*dWdx + dVdz.*dWdy));

% Render 3D isosurfaces
figure('Color', 'w');
Q_level = 0.1; % Adjust based on flow energy
p_iso = patch(isosurface(XI, YI, ZI, Q, Q_level));
p_iso.FaceColor = [0.9 0.65 0.1]; 
p_iso.EdgeColor = 'none';
axis equal; grid on;
view([-0.5 0.5 1.5]); 
camlight('headlight'); lighting gouraud;

xlim([-0.5, 3.0]);   % x-axis range [x_min, x_max]
ylim([-0.5, 0.5]);   % y-axis range [y_min, y_max]
zlim([0.0, 1.0]);    % z-axis range [z_min, z_max] (spanwise)

xlabel('x/c'); ylabel('y/c'); zlabel('z/c');
title('Coherent Structures (Q-Criterion)');

%% spanwise averages
% Average across the z-dimension (3rd index of meshgrid)
UQ_span_avg = mean(UQ, 3);
UQ_span_avg(UQ_span_avg > clim_val) = clim_val; UQ_span_avg(UQ_span_avg < -clim_val) = -clim_val;

VQ_span_avg = mean(VQ, 3);
VQ_span_avg(VQ_span_avg > clim_val) = clim_val; VQ_span_avg(VQ_span_avg < -clim_val) = -clim_val;

WQ_span_avg = mean(WQ, 3);
WQ_span_avg(WQ_span_avg > clim_val) = clim_val; WQ_span_avg(WQ_span_avg < -clim_val) = -clim_val;

vort_z_span_avg = mean(vort_z, 3);
vort_z_span_avg(vort_z_span_avg < -clim_val) = -clim_val;

%% Plot 2D slice of spanwise-averaged vorticity (omega_z)
figure('Color', 'w');
contourf(XI(:,:,1), YI(:,:,1), vort_z_span_avg, v_levels, 'LineStyle', 'none');
colormap(cmap); colorbar; axis equal;
xlim([-0.5 5]); ylim([-0.5 0.5])
title('Spanwise-Averaged Vorticity \langle\omega_z\rangle_z');

%% Plot 2D slice of spanwise-averaged velocity (u)
figure('Color', 'w');
contourf(XI(:,:,1), YI(:,:,1), UQ_span_avg, v_levels, 'LineStyle', 'none');
colormap(cmap); colorbar; axis equal;
title('Spanwise-Averaged Velocity \langle $u$ \rangle','Interpreter','LaTeX');

% Plot 2D slice of spanwise-averaged velocity (v)
figure('Color', 'w');
contourf(XI(:,:,1), YI(:,:,1), VQ_span_avg, v_levels, 'LineStyle', 'none');
colormap(cmap); colorbar; axis equal;
title('Spanwise-Averaged Velocity $v$');

% Plot 2D slice of spanwise-averaged velocity (w)
figure('Color', 'w');
contourf(XI(:,:,1), YI(:,:,1), WQ_span_avg, v_levels, 'LineStyle', 'none');
colormap(cmap); colorbar; axis equal;
title('Spanwise-Averaged Velocity $w$');