% strouhal based on aerodynamic coefficients on airfoil
clear; clc; 

%% --- 1. User-Defined Parameters ---
c = 1.0;        % chord length
U = 1.0;        % Free-stream velocity
% 1. Open the file
% 1. Open the file
filename = 'N0012_R4e3_A10_forces.txt';

% 2. Load Data 
% Read the file, skipping the header line
data = readmatrix(filename, 'NumHeaderLines', 1);

time = data(:, 1); % Time column 
Cd   = 2.*data(:, 2); % Drag force column 
Cl   = 2.*data(:, 3); % Lift force column

% --- 3. Pre-processing ---
% Find the steady-state portion
% For this script, we take the last 50% of the simulation time
start_idx = round(length(time) * 0.50);
t_steady  = time(start_idx:end);
Cl_steady = Cl(start_idx:end);
Cd_steady = Cd(start_idx:end);


% Calculate sampling frequency and interval
dt = mean(diff(t_steady));
Fs = 1/dt;
N  = length(Cl_steady);

% --- 4. Spectral Analysis (FFT) ---
% Remove mean (DC offset) to focus on oscillations
Cl_detrend = Cl_steady - mean(Cl_steady);

% Perform FFT
Y = fft(Cl_detrend);
P2 = abs(Y/N);            % Two-sided spectrum (both pos. & neg. frequencies
P1 = P2(1:floor(N/2)+1);  % Single-sided spectrum (0Hz to Nyquist freq (Fs/2))
P1(2:end-1) = 2*P1(2:end-1); % multiply by two for energy conservation
f = Fs*(0:(N/2))/N;       % Frequency vector 

% --- 5. Find Dominant Frequency & Strouhal Number ---
[max_amplitude, max_idx] = max(P1);
f_shedding = f(max_idx);
St = (f_shedding * c) / U;

% --- 6. Results & Visualization ---
fprintf('Shedding Frequency: %.4f Hz\n', f_shedding);
fprintf('Calculated Strouhal Number: %.4f\n', St);

figure;
subplot(2,1,1);
plot(t_steady, Cl_steady);
title('Statistically Steady-State Lift Coefficient History','FontSize',16);
xlabel('Time $[s]$','FontSize',16); ylabel('$C_L$','Rotation',0,'FontSize',16); grid on; 
ylim([0 2*max(Cl_steady)]);

subplot(2,1,2);
plot(f, P1);
title('Lift Coefficient Spectrum','FontSize',16);
xlabel('Strouhal Number','FontSize',16); ylabel('$\Delta C_L$','FontSize',16,'Rotation',0);
xlim([0, f_shedding*5]); % zoom in on relevant frequencies
grid on; 

