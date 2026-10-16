% Cd and Cl plotting for second half of simulation time
% Cd and Cl averages over that time
clear; clc; 

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');

% load data
data = readmatrix('N0012_R4e3_A10_forces.txt', 'NumHeaderLines', 1);
time = data(:,1);
Fx = data(:,2); 
Fy = data(:,3);

% define Simulation Constants
alpha_deg = 10.0;            % Angle of attack
alpha = deg2rad(alpha_deg);

Cd = 2*Fx/(1*1^2); 
Cl = 2*Fy/(1*1^2);

% only keep steady state (second half of simulation)
t_p = time(round(length(time)/2) : end);
Cl_p = Cl(round(length(time)/2) : end);
Cd_p = Cd(round(length(time)/2) : end);

% Mean values
fprintf('Mean Cl: %.3f\n', mean(Cl_p));
fprintf('Mean Cd: %.3f\n', mean(Cd_p));
fprintf('CL/CD ratio: %.3f\n', mean(Cl_p)/mean(Cd_p))

% amplitudes
amp_CL = mean([max(Cl_p)-mean(Cl_p), mean(Cl_p)-min(Cl_p)]);
amp_CD = mean([max(Cd_p)-mean(Cd_p), mean(Cd_p)-min(Cd_p)]);
disp(amp_CL); disp(amp_CD);

%% Plotting
figure('Color', 'w');
plot(t_p, Cl_p, 'r', 'LineWidth', 1.5, 'DisplayName', '$C_L$'); hold on;
yline(mean(Cl_p), 'g', 'LineWidth', 1.5, 'DisplayName', '$\overline{C_L}$'); 
yline(mean(Cd_p), 'k', 'LineWidth', 1.5, 'DisplayName', '$\overline{C_D}$');
plot(t_p, Cd_p, 'b', 'LineWidth', 1.5, 'DisplayName', '$C_D$');
grid on;
xlabel('Non-dimensional Time'); ylabel('Coefficient Value');
title(['NACA 4412 Aerodynamic Coefficients at $\alpha = $', num2str(alpha_deg), '$^\circ$']);
legend('Location', 'best'); xlim([t_p(1) t_p(end)]); ylim([0 1.5*max(Cl_p)]);
%yline(0,'LineWidth',1.5,'DisplayName', '0');

