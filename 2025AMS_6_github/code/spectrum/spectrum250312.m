% Clear environment
clear; clc; close all;

% Set parameters
E0 = 1;                        % Amplitude
w0 = 2.9767e15;                % Central angular frequency (Hz)
T = 1e-15;                     % Time width (s)
Omega = 2/T;                   % Corresponding frequency bandwidth (Hz)

% Define frequency range
w = linspace(w0 - 5*Omega, w0 + 5*Omega, 1000); 

% Compute the electric field in the frequency domain
E_w = E0 * exp(-((w - w0).^2) / Omega^2);

% Compute the electric field in the time domain using Fourier Transform
t = linspace(-5*T, 5*T, 1000); % Time range
E_t = ifftshift(ifft(E_w));    % Inverse Fourier Transform

% Plot frequency domain electric field
figure;
subplot(2,1,1);
plot(w, abs(E_w), 'b', 'LineWidth', 2);
xlabel('Frequency \omega (Hz)');
ylabel('|{E}(\omega)|');
title('Electric Field Distribution in Frequency Domain');
grid on;

% Plot time domain electric field
subplot(2,1,2);
plot(t, abs(E_t), 'r', 'LineWidth', 2);
xlabel('Time t (s)');
ylabel('|E(t)|');
title('Electric Field Distribution in Time Domain');
grid on;
%%

% Clear environment
clear; clc; close all;

% Set parameters
E0 = 1;                        % Amplitude
w0 = 2.9767e15;                % Central angular frequency (Hz)
T = 10e-15;                     % Time width (s)
Omega = 2/T;                   % Corresponding frequency bandwidth (Hz)

% Define frequency range
w = linspace(w0 - 5*Omega, w0 + 5*Omega, 1000); 

% Compute the electric field in the frequency domain
E_w = E0 * exp(-((w - w0).^2) / Omega^2);

% Compute the electric field in the time domain using Fourier Transform
t = linspace(-5*T, 5*T, 1000); % Time range
E_t = ifftshift(ifft(E_w));    % Inverse Fourier Transform

% Plot frequency domain electric field
figure(1);
plot(w, abs(E_w), 'b', 'LineWidth', 2);
xlabel('Frequency \omega (Hz)');
ylabel('|{E}(\omega)|');
title('Electric Field Distribution in Frequency Domain');
grid on;