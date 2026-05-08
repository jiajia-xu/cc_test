%% 1 Ex and Ey （single L）

% Constants
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.51078;
L = 1e-5;

% Dispersion parameters
dn_dlambda = -0.0198418e6;        % /m
d2n_dlambda2 = 0.0492482e12;      % /m^2
v_phi = c / n;        % Phase velocity [m/s]
k = omega0 / v_phi;   % Wavevector [rad/m]
k1 = (1/c) * (n + lambda0 * dn_dlambda);      % group velocity
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2; % group velocity delay

% Pulse parameters
t0 = 2e-14;      %T = 20fs Figure1(b)
Omega = 2 / t0;  %below eq.11
E0tilde = 1;

% Time array
t3 = linspace(-3*t0, 7*t0, 40000);  % Larger window

% A0 and input field
A0 = E0tilde * Omega * sqrt(pi); % below eq.11
Ex = A0 * exp(-t3.^2 / t0^2) .* exp(-1i * omega0 * t3) .* exp(1i * k0 * L);
Ex = real(Ex);

% Output field (analytic after GDD)
A = 1 / Omega^2 - 1i * k2 * L / 2;
B = k1 * L - t3;
Ey = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3) .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4*A));
Ey = real(Ey); 

% Plot
figure;
plot(t3 * 1e15, Ex, 'LineWidth', 1.2); hold on;
plot(t3 * 1e15, Ey, 'LineWidth', 1.2);
xlabel('Time (fs)');
ylabel('Electric Field');
legend('Ex (input)', 'Ey (after BK7)');
title('Ultrashort Pulse Propagation Through BK7 (Analytic GDD Only)');
grid on;

%% 2Good Ex and Ey （3 L）
% Constants
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.5108;

% Dispersion parameters
dn_dlambda = -0.0198e6;        % /m
d2n_dlambda2 = 0.0491e12;      % /m^2

% Pulse parameters
t0 = 2e-14; %T = 20fs
Omega = 2 / t0;
E0tilde = 1;

% Time array
t3 = linspace(-3*t0, 7*t0, 40000);

% A0 and input field
A0 = E0tilde * Omega * sqrt(pi);
Ex = A0 * exp(-t3.^2 / t0^2) .* exp(-1i * omega0 * t3);
Ex = real(Ex);

% custom_colors = [
%     0.890, 0.494, 0.478;  % 红色
%     0.467, 0.737, 0.494;  % 绿色
%     %0.933, 0.866, 0.510;  % 黄色
%     0.690, 0.478, 0.631; %紫色
% ];
custom_colors = [
    0.0, 0.9, 0.9  % Cyan
    0, 0, 1;   % Blue
    1, 0, 0;   % Red
];

% Thickness values for L
L_values = [0.5e-5, 1e-5, 1.5e-5];
%colors = ['r', 'g', 'b'];  % Color for each L

% Plot input field
figure('Position', [10 10 900 400]);
plot(t3 * 1e15, Ex, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5); hold on;%[0.302 0.737 0.890]

% Loop over each L
for i = 1:length(L_values)
    L = L_values(i);

    v_phi = c / n;
    k = omega0 / v_phi;
    k1 = (1/c) * (n + lambda0 * dn_dlambda);
    k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

    A = 1 / Omega^2 - 1i * k2 * L / 2;
    B = k1 * L - t3;
    Ey = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3) .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4*A));
    Ey = real(Ey);

    plot(t3 * 1e15, Ey, 'Color', custom_colors(i, :), 'LineWidth', 1.2);
end

grid on;
set(gca, 'FontSize', 14);        % 设置坐标轴字体 'FontName', 'Times New Roman',
%set(gcf, 'Color', 'w');                                         % 白色背景，导出图更干净
xlabel('$t$ [fs]', 'Interpreter', 'latex', 'FontSize', 16, 'FontName', 'Times New Roman');
ylabel('$E(t)$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName', 'Times New Roman');
%title('Ultrashort Pulse Propagation Through BK7', 'Interpreter', 'latex', 'FontSize', 16, 'FontName', 'Times New Roman'); % (Analytic GDD Only)
legend('$E_x$', ...
       '$E_y$, $L = 5\,\mu\mathrm{m}$', ...
       '$E_y$, $L = 10\,\mu\mathrm{m}$', ...
       '$E_y$, $L = 15\,\mu\mathrm{m}$', ...
       'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');

exportgraphics(gcf, 'BK7_Ex_Ey_LLL.png', 'Resolution', 300)

%% 2 Good Ex and Ey（3L）- 250728Normalized

% Constants
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.51078;

% Dispersion parameters
dn_dlambda = -0.0198418e6;        % /m
d2n_dlambda2 = 0.0492482e12;      % /m^2

% Pulse parameters
t0 = 2e-14; %T = 20fs
Omega = 2 / t0;
%E0tilde = 1;

% Time array
t3 = linspace(-3*t0, 7*t0, 40000);

% A0 and input field
%A0 = E0tilde * Omega * sqrt(pi);
A0 = 1/2;
E0tilde = A0/(Omega * sqrt(pi));
Ex = A0 * exp(-t3.^2 / t0^2) .* exp(-1i * omega0 * t3);
Ex = 2*real(Ex);

% custom_colors = [
%     0.890, 0.494, 0.478;  % 红色
%     0.467, 0.737, 0.494;  % 绿色
%    %0.933, 0.866, 0.510;  % 黄色
%     0.690, 0.478, 0.631;  %紫色
% ];
custom_colors = [
    0.0, 0.9, 0.9  % Cyan
    0, 0, 1;   % Blue
    1, 0, 0;   % Red
];

% Thickness values for L
L_values = [0.5e-5, 1e-5, 1.5e-5];
colors = ['r', 'g', 'b'];  % Color for each L

% Plot input field
figure('Position', [10 10 900 400]);
plot(t3 * 1e15, Ex, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5); hold on; %[0.302 0.737 0.890]

% Loop over each L
for i = 1:length(L_values)
    L = L_values(i);

    v_phi = c / n;
    k = omega0 / v_phi;
    k1 = (1/c) * (n + lambda0 * dn_dlambda);
    k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

    A = 1 / Omega^2 - 1i * k2 * L / 2;
    B = k1 * L - t3;
    Ey = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3) .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4*A));
    Ey = 2 * real(Ey);

    plot(t3 * 1e15, Ey, 'Color', custom_colors(i, :), 'LineWidth', 1.2);
end


ax = gca;
grid on
ax.Color = 'w';                % 背景白
ax.GridColor = [0.6 0.6 0.6];  % 深灰网格
ax.GridAlpha = 0.6;            % 半透明
ax.LineWidth = 1;              % 加粗坐标轴

set(gca, 'FontSize', 15);        % 设置坐标轴字体 'FontName', 'Times New Roman',
%set(gcf, 'Color', 'w');                                         % 白色背景，导出图更干净
xlabel('time [fs]', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('$E(t)$ [a.u.]', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times New Roman');
%title('Ultrashort Pulse Propagation Through BK7', 'Interpreter', 'latex', 'FontSize', 16, 'FontName', 'Times New Roman'); % (Analytic GDD Only)
legend('$E_x$', ...
       '$E_y$, $L = 5\,\mu\mathrm{m}$', ...
       '$E_y$, $L = 10\,\mu\mathrm{m}$', ...
       '$E_y$, $L = 15\,\mu\mathrm{m}$', ...
       'Interpreter', 'latex', 'FontSize', 16, 'FontName', 'Times New Roman');

exportgraphics(gcf, 'BK7_Ex_Ey_LLL_Normalized.png', 'Resolution', 300)

%% 3Test
% Constants
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.51078;

% Dispersion
dn_dlambda = -0.0198418e6;        % /m
d2n_dlambda2 = 0.0492482e12;      % /m^2

% Pulse parameters
t0 = 2e-14;
Omega = 2 / t0;
A0 = 1/2;
E0tilde = A0 / (Omega * sqrt(pi));

% Time values of interest
t_vals = [7e-15, 20e-15, 92e-15];  % in seconds

% Compute Ex at 7 fs and 20 fs
Ex_vals = zeros(1, 2);
for i = 1:2
    t = t_vals(i);
    Ex_complex = A0 * exp(-t.^2 / t0^2) .* exp(-1i * omega0 * t);
    Ex_vals(i) = 2 * real(Ex_complex);
end

% Dispersion parameters
L = 5e-6;  % 5 µm
v_phi = c / n;
k = omega0 / v_phi;
%k1 = 5.089e-9;
k1 = (1/c) * (n - lambda0 * dn_dlambda);
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

% Compute Ey at 30 fs
t = t_vals(3);
A = 1 / Omega^2 - 1i * k2 * L / 2;
AA = sqrt(A);
B = k1 * L - t;
Ey_complex = E0tilde * exp(1i * k * L) * exp(-1i * omega0 * t) * sqrt(pi / A) * exp(-B^2 / (4*A));
Ey_val = 2 * real(Ey_complex);

% Display results
fprintf('Ex at 7 fs  = %.6f\n', Ex_vals(1));
fprintf('Ex at 20 fs = %.6f\n', Ex_vals(2));
fprintf('Ey at 92 fs (L = 15 µm) = %.6f\n', Ey_val);

%% 3Test fs um
clear all

c = 0.3; %um/fs
lambda0 = 0.8; %um
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.51078;

dn_dlambda = -0.0198418; %/um
d2n_dlambda2 = 0.0492482; %/um^2

% Pulse parameters
t0 = 20; %fs
Omega = 2 / t0;
A0 = 1/2;
E0tilde = A0 / (Omega * sqrt(pi));

% Time values of interest
t_vals = [7, 20, 90]; 

Ex_vals = zeros(1, 2);
for i = 1:2
    t = t_vals(i);
    Ex_complex = A0 * exp(-t.^2 / t0^2) .* exp(-1i * omega0 * t);
    Ex_vals(i) = 2 * real(Ex_complex);
end


%fprintf('Ex at 7 fs  = %.6f\n', Ex_vals(1));
%fprintf('Ex at 20 fs = %.6f\n', Ex_vals(2));


L = 15; %um
v_phi = c / n;
k = omega0 / v_phi;%

k1 = (1/c) * (n - lambda0 * dn_dlambda);
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

t = t_vals(3);
A = 1 / Omega^2 - 1i * k2 * L / 2;
B = k1 * L - t;
Ey_complex = E0tilde * exp(1i * k * L) * exp(-1i * omega0 * t) * sqrt(pi / A) * exp(-B^2 / (4*A));
Ey_val = 2 * real(Ey_complex);


fprintf('Ey at 30 fs (L = 5 µm) = %.6f\n', Ey_val);

%% 4  good Ex and Ey（1L）- 250901Normalized

% Constants
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.51078;

% Dispersion parameters
dn_dlambda = -0.0198418e6;        % /m
d2n_dlambda2 = 0.0492482e12;      % /m^2


% Pulse parameters
t0 = 2e-14; %T = 20fs
Omega = 2 / t0;
%E0tilde = 1;

% Time array
t3 = linspace(1.615*10^3*t0, 1.64*10^3*t0, 10000);

% A0 and input field
%A0 = E0tilde * Omega * sqrt(pi);
L_values = 6.4e-3;
t_delay_x = 32568*10^(-15);
%t_delay_x =0;
A0 = 1/2;
E0tilde = A0/(Omega * sqrt(pi));
Ex = A0 * exp(-(t3-t_delay_x).^2 / t0^2) .* exp(-1i * omega0 * t3) ;%.* exp(1i * k * L) 
Ex = 2*real(Ex);

% custom_colors = [
%     0.890, 0.494, 0.478;  % 红色
%     0.467, 0.737, 0.494;  % 绿色
%    %0.933, 0.866, 0.510;  % 黄色
%     0.690, 0.478, 0.631;  %紫色
% ];
custom_colors = [
    0.0, 0.9, 0.9  % Cyan
    0, 0, 1;   % Blue
    1, 0, 0;   % Red
];

% Thickness values for L

%L_values = [0.5e-5, 1e-5, 1.5e-5];
%colors = ['r', 'g', 'b'];  % Color for each L

% Plot input field
figure('Position', [10 10 900 400]);
plot(t3 * 1e15, Ex, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5); hold on;%[0.302 0.737 0.890]

% Loop over each L
for i = 1:length(L_values)
    L = L_values(i);

    v_phi = c / n;
    k = omega0 / v_phi;
    k1 = (1/c) * (n - lambda0 * dn_dlambda);
    k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

    A = 1 / Omega^2 - 1i * k2 * L / 2;
    B = k1 * L - t3;
    Ey = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3) .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4*A));
    Ey = 2 * real(Ey);

   plot(t3 * 1e15, Ey, 'Color', custom_colors(i, :), 'LineWidth', 1.2);
end


grid on;
set(gca, 'FontSize', 14);        % 设置坐标轴字体 'FontName', 'Times New Roman',
%set(gcf, 'Color', 'w');                                         % 白色背景，导出图更干净
xlabel('Time [fs]', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('$E(t)$', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times New Roman');
%title('Ultrashort Pulse Propagation Through BK7', 'Interpreter', 'latex', 'FontSize', 16, 'FontName', 'Times New Roman'); % (Analytic GDD Only)
legend('$E_x$', ...
       '$E_y$, $L = 6.4\,\mathrm{mm}$','Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');

%exportgraphics(gcf, 'BK7_Ex_Ey_LLL_Normalized.png', 'Resolution', 300)



%% 4  good Ex and Ey（1L）- 260312 Normalized 

% Constants
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.51078;

% Dispersion parameters
dn_dlambda = -0.0198418e6;        % /m
d2n_dlambda2 = 0.0492482e12;      % /m^2


% Pulse parameters
t0 = 2e-14; %T = 20fs
Omega = 2 / t0;
%E0tilde = 1;

% Time array
t3 = linspace(1.615*10^3*t0, 1.64*10^3*t0, 10000);

% A0 and input field
%A0 = E0tilde * Omega * sqrt(pi);
L_values = 6.4e-3;
t_delay_x = 32568*10^(-15);
%t_delay_x =0;
A0 = 1/2;
E0tilde = A0/(Omega * sqrt(pi));
Ex = A0 * exp(-(t3-t_delay_x).^2 / t0^2) .* exp(-1i * omega0 * t3) ;%.* exp(1i * k * L) 
Ex = 2*real(Ex);

% custom_colors = [
%     0.890, 0.494, 0.478;  % 红色
%     0.467, 0.737, 0.494;  % 绿色
%    %0.933, 0.866, 0.510;  % 黄色
%     0.690, 0.478, 0.631;  %紫色
% ];
custom_colors = [
    0.0, 0.9, 0.9  % Cyan
    0, 0, 1;   % Blue
    1, 0, 0;   % Red
];

% Thickness values for L

%L_values = [0.5e-5, 1e-5, 1.5e-5];
%colors = ['r', 'g', 'b'];  % Color for each L

% Plot input field
figure('Position', [10 10 900 400]);
plot(t3 * 1e15, Ex, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5); hold on;%[0.302 0.737 0.890]

% Loop over each L
for i = 1:length(L_values)
    L = L_values(i);

    v_phi = c / n;
    k = omega0 / v_phi;
    k1 = (1/c) * (n - lambda0 * dn_dlambda);
    k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

    A = 1 / Omega^2 - 1i * k2 * L / 2;
    B = k1 * L - t3;
    Ey = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3) .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4*A));
    Ey = 2 * real(Ey);

   plot(t3 * 1e15, Ey, 'Color', custom_colors(i, :), 'LineWidth', 1.2);
end


grid on;
set(gca, 'FontSize', 14);        % 设置坐标轴字体 'FontName', 'Times New Roman',
%set(gcf, 'Color', 'w');                                         % 白色背景，导出图更干净
xlabel('Time [fs]', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('$E(t)$', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times New Roman');
%title('Ultrashort Pulse Propagation Through BK7', 'Interpreter', 'latex', 'FontSize', 16, 'FontName', 'Times New Roman'); % (Analytic GDD Only)
legend('$E_x$', ...
       '$E_y$, $L = 6.4\,\mathrm{mm}$','Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');

%exportgraphics(gcf, 'BK7_Ex_Ey_LLL_Normalized.png', 'Resolution', 300)


%% Ex and Ey（1L）- 250904
% Constants
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.51078;

% Dispersion parameters
dn_dlambda = -0.0198418e6;        % /m
d2n_dlambda2 = 0.0492482e12;      % /m^2

% Pulse parameters
t0 = 2e-14; %T = 20fs
Omega = 2 / t0;

% Time array
t3 = linspace(3.24/2*10^3*t0, 3.275/2*10^3*t0, 5000);

% Propagation parameters
L_values = 6.4e-3;  % propagation distance

% Initial amplitude
A0 = 1/2;
E0tilde = A0/(Omega * sqrt(pi));

% Colors
custom_colors = [
    0.0, 0.9, 0.9;   % Cyan
    0, 0, 1;         % Blue
    1, 0, 0;         % Red
];

% Precompute dispersion terms
v_phi = c / n;
k = omega0 / v_phi;
k1 = (1/c) * (n - lambda0 * dn_dlambda);
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

% Delay list
t_delay_list = linspace(32490, 32650, 161) * 1e-15;

% ---- Create video ----
v = VideoWriter('pulse_propagation.mp4','MPEG-4');
v.FrameRate = 20;   % frame rate (fps)
open(v);

figure('Position', [10 10 900 400]);

for j = 1:length(t_delay_list)
    t_delay_x = t_delay_list(j);

    % Input field
    Ex = A0 * exp(-(t3-t_delay_x).^2 / t0^2) .* exp(-1i * omega0 * t3);
    Ex = 2*real(Ex);

    % Output field
    L = L_values;
    A = 1 / Omega^2 - 1i * k2 * L / 2;
    B = k1 * L - t3;
    Ey = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3) ...
         .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4*A));
    Ey = 2 * real(Ey);

    % ---- Plot ----
    clf;
    plot(t3*1e15, Ex, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5); hold on;
    plot(t3*1e15, Ey, 'Color', custom_colors(2,:), 'LineWidth', 1.2);
    grid on;
    set(gca,'FontSize',14);
    xlabel('$t$ [fs]', 'Interpreter','latex','FontSize',16);
    ylabel('$E(t)$', 'Interpreter','latex','FontSize',16);
    %legend('$E_x$', sprintf('$E_y$, L = %.1f mm', L*1e3), ...'Interpreter','latex','FontSize',14,'Location','best');
    legend('$E_x$', '$E_y$', 'Interpreter','latex','FontSize',14,'Location','best');
    title(sprintf('$\\Delta t = %.1f\\,\\mathrm{fs}$', t_delay_x*1e15), ...
          'Interpreter','latex','FontSize',16);

    % Save frame
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
disp('Video saved: pulse_propagation.mp4');


%% 5 checking 260313- instantaneous Stokes parameters
clear; clc; close all;

% Constants
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
n = 1.51078;

% Dispersion parameters (BK7 at 800 nm)
dn_dlambda   = -0.0198418e6;    % /m
d2n_dlambda2 =  0.0492482e12;   % /m^2

% Pulse parameters
t0 = 20e-15;                    % 20 fs
Omega = 2 / t0;

% Time array
t = linspace(-40e-15, 70e-15, 40000);   % [-40,70] fs, like your example

% Input field amplitude
A0 = 1/2;
E0tilde = A0/(Omega * sqrt(pi));

% Input x-component (complex analytic field)
Ex_c = A0 * exp(-t.^2 / t0^2) .* exp(-1i * omega0 * t);

% Only one length: L = 4 um
L = 4e-6;

% Dispersion quantities
v_phi = c / n;
k  = omega0 / v_phi;
k1 = (1/c) * (n + lambda0 * dn_dlambda);
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

A = 1 / Omega^2 - 1i * k2 * L / 2;
B = k1 * L - t;

% Output y-component (complex analytic field)
Ey_c = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t) ...
     .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4*A));

% Instantaneous Stokes parameters
S0 = abs(Ex_c).^2 + abs(Ey_c).^2;

% avoid division by zero in the far tails
eps0 = 1e-14 * max(S0);
S0_safe = S0;
S0_safe(S0_safe < eps0) = NaN;

s1 = (abs(Ex_c).^2 - abs(Ey_c).^2) ./ S0_safe;
s2 = 2 * real(Ex_c .* conj(Ey_c)) ./ S0_safe;
s3 = 2 * imag(Ex_c .* conj(Ey_c)) ./ S0_safe;
% 如果你发现 s3 的正负和参考图相反，改成：
% s3 = -2 * imag(Ex_c .* conj(Ey_c)) ./ S0_safe;

% Normalized intensity for plotting
Int = S0 / max(S0);

% Plot
figure('Position',[100 100 820 500]); hold on;

plot(t*1e15, Int, 'k--', 'LineWidth', 2.0);
plot(t*1e15, s1,  'b',   'LineWidth', 2.0);
plot(t*1e15, s2,  'Color', [1.0 0.5 0.0], 'LineWidth', 2.0);
plot(t*1e15, s3,  'g',   'LineWidth', 2.0);

grid on;
box on;
set(gca, 'FontSize', 14, 'LineWidth', 1.0);
xlabel('time [fs]', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Stokes parameters', 'Interpreter', 'latex', 'FontSize', 16);
xlim([-40 70]);
ylim([-1.05 1.25]);

legend('Int', '$s_1$', '$s_2$', '$s_3$', ...
    'Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');

title('Evolution of the three instantaneous Stokes parameters ($L=4\,\mu$m)', ...
    'Interpreter', 'latex', 'FontSize', 15);

exportgraphics(gcf, 'BK7_Stokes_L4um.png', 'Resolution', 300);


%% 6 checking260313
clear; clc; close all;

% Constants
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
n = 1.51078;

% Dispersion parameters
dn_dlambda   = -0.0198418e6;    % /m
d2n_dlambda2 =  0.0492482e12;   % /m^2

% Pulse parameters
t0 = 20e-15;                    % 20 fs
Omega = 2 / t0;

% Time array
t = linspace(-40e-15, 70e-15, 40000);

% Input amplitude
A0 = 1/2;
E0tilde = A0/(Omega * sqrt(pi));

% Thickness values
L_values = [2e-6, 4e-6, 6e-6, 8e-6];

% Dispersion quantities
v_phi = c / n;
k  = omega0 / v_phi;
k1 = (1/c) * (n + lambda0 * dn_dlambda);
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

% Figure
figure('Position',[100 100 820 470]); hold on;

for i = 1:length(L_values)
    L = L_values(i);

    A = 1/Omega^2 - 1i*k2*L/2;
    B = k1*L - t;

    % complex analytic output field
    Ey_c = E0tilde .* exp(1i*k*L) .* exp(-1i*omega0*t) ...
         .* sqrt(pi./A) .* exp(-B.^2./(4*A));

    Ex_c = A0 * exp(-t.^2 / t0^2) .* exp(-1i * omega0 * t);
    I = abs(Ex_c).^2 + abs(Ey_c).^2;
   
    % normalize each curve by its own maximum
    %I = I / max(I);

    plot(t*1e15, I, 'LineWidth', 1.8);
end

% Axis style
grid on;
box on;
ax = gca;
ax.LineWidth = 1.0;
ax.FontSize = 13;
ax.Color = [0.97 0.97 0.97];
ax.GridColor = [0.65 0.65 0.65];
ax.GridAlpha = 0.55;

xlabel('time [fs]', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Intensity [a.u.]', 'Interpreter', 'latex', 'FontSize', 16);

xlim([-40 70]);
ylim([0 0.5]);

legend({'$2\,\mu\mathrm{m}$', '$4\,\mu\mathrm{m}$', '$6\,\mu\mathrm{m}$', '$8\,\mu\mathrm{m}$'}, ...
       'Interpreter', 'latex', 'FontSize', 13, 'Location', 'northeast');

exportgraphics(gcf, 'BK7_Intensity_Comparison_4L.png', 'Resolution', 300);