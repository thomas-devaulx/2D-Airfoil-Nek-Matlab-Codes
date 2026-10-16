%% Nek5000 Post-Processing: Combined Vorticity + Mesh Overlay
clear; clc; close all;
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 

%% Load Data

[data] = readnek('N0012_R4e3_A10_dd0.f00334');

X = data(:,:,1); xg = X(:);
Y = data(:,:,2); yg = Y(:);
u = data(:,:,3); ug = u(:);
v = data(:,:,4); vg = v(:);

[n_elem, pts_per_elem] = size(X);
N = sqrt(pts_per_elem); 

% Create Velocity Interpolants
F_u = scatteredInterpolant(xg, yg, ug, 'linear', 'none');
F_v = scatteredInterpolant(xg, yg, vg, 'linear', 'none');

% formula: 4*N - 2
pts_per_outline = 4*N - 2; 
outline_x = zeros(pts_per_outline, n_elem);
outline_y = zeros(pts_per_outline, n_elem);

for e = 1:n_elem
    x_e = reshape(X(e,:), [N, N]);
    y_e = reshape(Y(e,:), [N, N]);
    
    % Trace 4 edges of the element, returning to start, and separate with NaN
    ox = [x_e(:,1)', x_e(N, 2:N), x_e(N-1:-1:1, N)', x_e(1, N-1:-1:1), NaN];
    oy = [y_e(:,1)', y_e(N, 2:N), y_e(N-1:-1:1, N)', y_e(1, N-1:-1:1), NaN];
    
    outline_x(:, e) = ox';
    outline_y(:, e) = oy';
end

mesh_x = outline_x(:);
mesh_y = outline_y(:);

%% Visualization Parameters
clim_val = 20; clim_val_p = 5; % limits of the vorticity field colorbar
v_levels = linspace(-clim_val, clim_val, 120); % number of levels, (., ., 60)
cmap = b2r(-clim_val, clim_val); 

% Figure 1: full zoomed out wake

nx = 3200; ny = 1000;
xi = linspace(-1.0, 40.0, nx); yi = linspace(-2.0, 2.0, ny); % provides values over which the vorticity is calculated
[XI, YI] = meshgrid(xi, yi);
UQ = F_u(XI, YI); VQ = F_v(XI, YI);
[~, dUdy] = gradient(UQ, xi, yi); [dVdx, ~] = gradient(VQ, xi, yi);
vortQ = dVdx - dUdy; vortQ(vortQ > clim_val) = clim_val; vortQ(vortQ < -clim_val) = -clim_val;

% Plotting
figure('Color', 'w', 'Name', 'Wake View', 'Position', [100, 25, 1200, 400]);
contourf(XI, YI, vortQ, v_levels, 'LineStyle', 'none'); hold on;
colormap(cmap); clim([-clim_val_p, clim_val_p]); colorbar;
axis equal; xlim([-1, 40]); ylim([-2, 2]); % domain over which the plot will initially be focused on
title('Far-Wake');
xlabel('x/c'); ylabel('y/c');


%% Vortex structure label and pinpoint

% get max and min values of the calculated vorticity
[max_vort_values, max_id_y_all] = max(vortQ); 
[min_vort_values, min_id_y_all] = min(vortQ); 

peak_id_x = []; trough_id_x = [];

% find the associated id of the max values and write to a new array is the max is a max of maxes
for i = 2:length(max_vort_values)-1
    if max_vort_values(i) > max_vort_values(i-1) && max_vort_values(i) > max_vort_values(i+1)
        peak_id_x = [peak_id_x; i];
    end
end
for i = 2:length(min_vort_values)-1
    if min_vort_values(i) < min_vort_values(i-1) && min_vort_values(i) < min_vort_values(i+1)
        trough_id_x = [trough_id_x; i];
    end
end

% Peaks of magnitude > 3
keep_peaks = abs(max_vort_values(peak_id_x)) > 3; 
peak_id_x = peak_id_x(keep_peaks);                  
maximum_vort = max_vort_values(peak_id_x);          

% Troughs of magnitude > 3
keep_troughs = abs(min_vort_values(trough_id_x)) > 3; 
trough_id_x = trough_id_x(keep_troughs);                
minimum_vort = min_vort_values(trough_id_x);            

% finds corresponding coordinate values from the x id
peak_x_coord = peak_id_x * 41/nx - 1; 
trough_x_coord = trough_id_x * 41/nx - 1;

% ensuring only keep the vortices after a certain distance
pos_peaks_mask   = peak_x_coord >= 5;
pos_troughs_mask = trough_x_coord >= 5;

peak_id_x    = peak_id_x(pos_peaks_mask);
maximum_vort = maximum_vort(pos_peaks_mask);
peak_x_coord = peak_x_coord(pos_peaks_mask);
trough_id_x    = trough_id_x(pos_troughs_mask);
minimum_vort   = minimum_vort(pos_troughs_mask);
trough_x_coord = trough_x_coord(pos_troughs_mask);

% get the y values for the peaks and troughs
peak_id_y   = max_id_y_all(peak_id_x);
trough_id_y = min_id_y_all(trough_id_x);

% convert to coordinates
peak_y_coord   = peak_id_y * 4/ny - 2; 
trough_y_coord = trough_id_y * 4/ny - 2;

%% Figure 2: multi level vorticity contours from x = -1
figure('Position', [100, 500, 1200, 400]);
hold on;

% vorticity contour values to draw these contour lines, add remove or replace
contour_levels = [0.1 0.5, 1.0, 1.5, 2.0, 2.5]; 

% Plot positive contours in red (x from -1 to 20)
[~, h_pos] = contour(XI, YI, vortQ, contour_levels, 'r', 'LineWidth', 1.0);

% Plot negative contours in blue (x from -1 to 20)
[~, h_neg] = contour(XI, YI, vortQ, -contour_levels, 'b', 'LineWidth', 1.0);

title('Multi-Level Wake Vorticity Contours with Targeted Highlights');
xlabel('x/c'); ylabel('y/c');
axis equal; xlim([-1.0 40.0]); ylim([-2.0 2.0]); grid on;
set(gca, 'Color', 'w', 'Layer', 'top', 'GridColor', [0.8 0.8 0.8]); 

%% Figure 3: Decay rate of core values, not log plot
figure('Color', 'w', 'Name', 'Vortex_Decay_Rate', 'Position', [150, 150, 700, 450]);
hold on;

% plot positive core decay
plot(peak_x_coord, maximum_vort, 'o-', 'LineWidth', 2, ...
     'MarkerSize', 6, 'MarkerFaceColor', 'r', 'Color', 'r', ...
     'DisplayName', 'Positive Cores (\omega_{max})');

% plot negative core decay
plot(trough_x_coord, minimum_vort, 'o-', 'LineWidth', 2, ...
     'MarkerSize', 6, 'MarkerFaceColor', 'b', 'Color', 'b', ...
     'DisplayName', 'Negative Cores Magnitude |\omega_{min}|');

yline(0, 'LineWidth', 2);
title('Peak Vorticity Strength Spatial Decay','FontSize',16);
xlabel('Downstream Core Position (x/c)','FontSize',16);
ylabel('Peak Vorticity Value $\frac{(\omega_{peak})c}{U_\infty}$','FontSize',16);
xlim([1.2 40.0]);
ylim([-clim_val clim_val]);
grid on;
set(gca, 'GridColor', [0.8 0.8 0.8]);

%% Figure 4: decay rate of core values, loglog plot
figure('Color', 'w', 'Name', 'Log_Log_Vortex_Decay_Rate', 'Position', [175, 175, 700, 450]);
hold on;

% plot positive core decay 
plot(peak_x_coord, maximum_vort, 'o-', 'LineWidth', 2, ...
     'MarkerSize', 6, 'MarkerFaceColor', 'r', 'Color', 'r', ...
     'DisplayName', 'Positive Cores ($\omega_{max}$)');

% plot negative core decay
plot(trough_x_coord, abs(minimum_vort), 'o-', 'LineWidth', 2, ...
     'MarkerSize', 6, 'MarkerFaceColor', 'b', 'Color', 'b', ...
     'DisplayName', 'Negative Cores Magnitude $|\omega_{min}|$');

% Extract data for the fit
x_fit = peak_x_coord(1:end);
y_fit = maximum_vort(1:end);
x_fit_neg = trough_x_coord(1:end);
y_fit_neg = abs(minimum_vort(1:end));
    
% polyfit: log10(y) = m*log10(x) + b
p_pos = polyfit(log10(x_fit), log10(y_fit), 1);
decay_exponent = p_pos(1); % m
intercept = p_pos(2);      % b

p_neg = polyfit(log10(x_fit_neg), log10(y_fit_neg), 1);
decay_exponent_neg = p_neg(1); % m
intercept_neg = p_neg(2);      % b
    
% Coordinates for the line
x_ref = linspace(min(x_fit), max(x_fit), 100);
y_ref = 10.^(decay_exponent * log10(x_ref) + intercept);

x_ref_neg = linspace(min(x_fit_neg), max(x_fit_neg), 100);
y_ref_neg = 10.^(decay_exponent_neg * log10(x_ref_neg) + intercept_neg);
    
% Plot 
plot(x_ref, y_ref, 'k', 'LineWidth', 1.8, 'HandleVisibility', 'off');
plot(x_ref_neg, y_ref_neg, 'Color', [1 0.65 0], 'LineWidth', 1.8, 'HandleVisibility', 'off');

idx_pos = round(length(x_ref) * 0.80);
idx_neg = round(length(x_ref_neg) * 0.80);
text(x_ref(idx_pos), y_ref(idx_pos)-2, sprintf('$\\sim x^{%.2f}$', decay_exponent), ...
     'Interpreter', 'latex', 'FontSize', 14, 'Color', 'k', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
text(x_ref_neg(idx_neg), y_ref_neg(idx_neg)+2.5, sprintf('$\\sim x^{%.2f}$', decay_exponent_neg), ...
     'Interpreter', 'latex', 'FontSize', 14, 'Color', [1 0.65 0], ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

set(gca, 'XScale', 'log', 'YScale', 'log');

title('Log-Log Peak Vorticity Value Spatial Decay','FontSize',16);
xlabel('Downstream Position (x/c)','FontSize',16);
ylabel('$\frac{|\omega_{peak}|c}{U_\infty}$','FontSize',16,'Rotation',0);


xlim([5 40.0]); ylim([.1 30]);
grid on;
set(gca, 'GridColor', [0.8 0.8 0.8]);
legend('Location', 'best','FontSize',16);


%% Figure 6: b/a and a/c spacing
c = 1;
a_peak = zeros(length(peak_x_coord),1);
a_trough = zeros(length(trough_x_coord),1);

for i = 2:length(peak_x_coord)
    a_peak(i) = peak_x_coord(i) - peak_x_coord(i-1);
end
for i = 2:length(trough_x_coord)
    a_trough(i) = trough_x_coord(i) - trough_x_coord(i-1);
end

ac_spacing = figure('Position',[250 250 1000 300],'Name','long_spac'); 
hold on;
plot(peak_x_coord(2:end-1), a_peak(2:end-1),'+r','LineWidth',2);
plot(trough_x_coord(2:end-1), a_trough(2:end-1),'+b','LineWidth',2);

ylim([0 5]); xlim([0 40]);
%yline(0.852); % edit to match the strouhal number
grid on;
legend('Positive Vortices','Negative Vortices','Location','SouthWest','FontSize',16);
xlabel('x/c','FontSize',16); ylabel('a/c','FontSize',16,'Rotation',0); 
title('Longitudinal Separation between vortices of the same sign','FontSize',16);

%%
b_peak = zeros(length(peak_y_coord),1);
for i = 1:length(peak_y_coord)-6
    b_peak(i) = abs(abs(trough_y_coord(i)) - abs(peak_y_coord(i)));
end
b_trough = zeros(length(peak_y_coord),1);
for i = 1:length(peak_y_coord)-6
    b_trough(i) = abs(abs(peak_y_coord(i)) - abs(trough_y_coord(i)));
end


ba_peak = b_peak ./ a_peak(1:end);
ba_trough = b_trough(1:end) ./ a_trough(3:end);

ba_spacing = figure('Position',[250 350 1000 400],'Name','lat_spac');
hold on;
plot(peak_x_coord(1:end-1), ba_peak(1:end-1),'+r','LineWidth',2);
plot(trough_x_coord(1:end-3), ba_trough(1:end-1),'+b','LineWidth',2)

xlim([0 40]); ylim([0 0.5]); xticks(0:1:40); yticks(0:0.05:0.5);
yline(0.281,'LineWidth',1.5);
grid on;
legend('Positive Vortices','Negative Vortices','von Karman Ratio, 0.281','Location','NorthWest','FontSize',16);
xlabel('x/c','FontSize',16); ylabel('b/a','FontSize',16,'Rotation',0); 
title('Lateral Spacing between opposing sign vortices','FontSize',16);