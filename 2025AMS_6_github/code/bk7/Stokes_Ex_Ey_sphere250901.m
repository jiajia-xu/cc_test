%% 1
clc; clear all;
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

v_phi = c / n;
k = omega0 / v_phi;%
L = 6400;
t_delay_x =32568; t = linspace(3.246/2*10^3*t0, 3.266/2*10^3*t0, 201);
%t_delay_x = 32584; t = linspace(3.25/2*10^3*t0, 3.27/2*10^3*t0, 201);
%t_delay_x = 32594;

k1 = (1/c) * (n - lambda0 * dn_dlambda);
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;




Ex_complex = A0 * exp(-(t-t_delay_x).^2 / t0^2) .* exp(-1i * omega0 * t).* exp(1i * k0 * L);
%Ex = 2 * real(Ex_complex);

A = 1 / Omega^2 - 1i * k2 * L / 2;
B = k1 * L - t;
Ey_complex = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t) .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4*A));
%Ey= 2 * real(Ey_complex);

S0 = abs(Ex_complex).^2 + abs(Ey_complex).^2;


% figure;
% set(gcf, 'Position', [100, 100, 1200, 400]);  % 1:3 横向图像比例
% %plot(t * 1e15, S0_norm, 'm', 'LineWidth', 2);
% plot(t , S0, 'm', 'LineWidth', 2);
% xlabel('Time [fs]', 'FontSize', 16);
% ylabel('S_0(t)', 'FontSize', 16);
% title('Stokes parameter S_0(t) = |E_x|^2 + |E_y|^2', 'FontSize', 18);
% grid on;
% set(gca, 'FontSize', 14);

%% 2 s2
S1 = abs(Ex_complex).^2 - abs(Ey_complex).^2;
s1 = S1./S0;

% 归一化到最大绝对值（保持正负对称）
%S1_norm = S1 / max(abs(S1));

% 绘图
figure;

set(gcf, 'Position', [100, 100, 1200, 400]);  % 横向图像比例 1:3
plot(t , s1, 'c', 'LineWidth', 2);
ylim([-1,1]);
set(gca, 'FontSize', 14);
xlabel('$t$ [fs]', 'Interpreter', 'latex', 'FontSize', 24, 'FontName', 'Times New Roman');
ylabel('$s1(t)$', 'Interpreter', 'latex', 'FontSize', 24, 'FontName', 'Times New Roman');
title('Normalized Stokes parameter $s_1(t)$', 'FontSize', 18, 'Interpreter', 'latex');


grid on;

%exportgraphics(gcf, 's1_L15um.png', Resolution=300);
%exportgraphics(gcf, 'S1_L_10um.png', 'Resolution', 300);



  

%% 3Four lines in one plots
% 计算 Stokes 参数
    S0 = abs(Ex_complex).^2 + abs(Ey_complex).^2;
    S1 = abs(Ex_complex).^2 - abs(Ey_complex).^2;
    %S2 = 2 * real(Ex_complex .* conj(Ey_complex));
    %S3 = - 2 * imag(Ex_complex .* conj(Ey_complex));
    S2 = conj(Ex_complex).*Ey_complex + conj(Ey_complex) .*Ex_complex;
    S3 = 1i .* [conj(Ey_complex) .*Ex_complex - conj(Ex_complex).*Ey_complex];

    s0 = (S1.^2+S2.^2+S3.^2)./S0.^2;
    s1 = S1./S0;
    s2 = S2./S0;
    s3 = S3./S0;

% 画图
figure;
set(gcf, 'Position', [100, 100, 1200, 400]);

hold on;

plot(t , s0, 'k', 'LineWidth', 2, 'DisplayName', 'sum $s^2$');
plot(t , s1, 'r', 'LineWidth', 2, 'DisplayName', '$s_1$');
plot(t , s2, 'g', 'LineWidth', 2, 'DisplayName', '$s_2$');
plot(t , s3, 'b', 'LineWidth', 2, 'DisplayName', '$s_3$');

ylim([-1,1]);
set(gca, 'FontSize', 14);
xlabel('$t$ [fs]', 'Interpreter', 'latex', 'FontSize', 16, 'FontName', 'Times New Roman');
ylabel('$s(t)$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName', 'Times New Roman');
%title('Four Stokes Parameters', 'FontSize', 18, 'FontName', 'Times New Roman');
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);
grid on;

hold off;
