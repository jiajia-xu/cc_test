%% 3 - Plot S_0(t) = |E_x|^2 + |E_y|^2 for given L
% 使用 L = 1.0e-5，对应前面图一中输入 E_x 和图二中输出 E_y
% Constants
c = 3e8;              % Speed of light [m/s]
n = 1.453;            % Refractive index for fused silica
v_phi = c / n;        % Phase velocity [m/s]
omega_0 = 2.36e15;    % Central angular frequency [rad/s]
lambda = 0.8e-6;      % Wavelength [m]
L = 3.0e-5;           % Propagation distance [m]
t0 = 1e-14;           % Pulse width [s]
k = omega_0 / v_phi;  % Wavevector [rad/m]
k1 = 5.15e-9;         % Group delay [s/m]
k2 = 1.99e-23;        % Dispersion
Omega = 2e14;         % Spectral width [rad/s]
E0tilde = 1;
A0 = E0tilde * Omega * sqrt(pi); 

% 时间轴：我们统一选最长时间范围以便兼容
t3 = linspace(-1e-13, 4e-13, 10000);

% 补齐 E_x(t): 使用与图1一致的表达式
E_x = A0 * exp(-t3.^2 / t0^2) .* exp(-1i * omega_0 * t3);

% E_y(t) for L = 1.0e-5（与图2保持一致）
A = 1/Omega^2 - 1i * k2 * L / 2;
B = k1 * L - t3;
E_y = E0tilde .* exp(1i * k * L) .* exp(-1i * omega_0 * t3) .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4 * A));


% 计算 S0(t)
S0 = abs(E_x).^2 + abs(E_y).^2;

% 归一化
%S0_norm = S0 / max(S0);

% 绘图
figure(4);
set(gcf, 'Position', [100, 100, 1200, 400]);  % 1:3 横向图像比例
%plot(t3 * 1e15, S0_norm, 'm', 'LineWidth', 2);
plot(t3 * 1e15, S0, 'm', 'LineWidth', 2);
xlabel('Time [fs]', 'FontSize', 16);
ylabel('S_0(t)', 'FontSize', 16);
title('Stokes parameter S_0(t) = |E_x|^2 + |E_y|^2', 'FontSize', 18);
grid on;
set(gca, 'FontSize', 14);

% 保存图像
exportgraphics(gcf, 'S0_L_30um.png', 'Resolution', 300);


%% 5 Stokes 参数 S_1(t) = |E_x|^2 - |E_y|^2 --------

% 计算 S1(t)
S1 = abs(E_x).^2 - abs(E_y).^2;

% 归一化到最大绝对值（保持正负对称）
%S1_norm = S1 / max(abs(S1));

% 绘图
figure(5);
set(gcf, 'Position', [100, 100, 1200, 400]);  % 横向图像比例 1:3
plot(t3 * 1e15, S1, 'c', 'LineWidth', 2);
%plot(t3 * 1e15, S1_norm, 'c', 'LineWidth', 2);
xlabel('Time [fs]', 'FontSize', 16);
ylabel('S_1(t)', 'FontSize', 16);
title('Stokes parameter S_1(t) = |E_x|^2 - |E_y|^2', 'FontSize', 18);
grid on;
set(gca, 'FontSize', 14);

% 保存图像
exportgraphics(gcf, 'S1_L_30um.png', 'Resolution', 300);

%% 6  Stokes 参数 S_2(t) = 2 * Re(E_x * conj(E_y)) --------
%S2 = real(conj(E_x) .* E_y) + real(conj(E_y) .* E_x);
S2 = 2 * real(E_x .* conj(E_y));
S2_norm = S2 / max(abs(S2));

figure(6);
set(gcf, 'Position', [100, 100, 1200, 400]);
%plot(t3 * 1e15, S2_norm, 'g', 'LineWidth', 2);
plot(t3 * 1e15, S2, 'g', 'LineWidth', 2);
xlabel('Time [fs]', 'FontSize', 16);
ylabel('S_2(t)', 'FontSize', 16);
title('Stokes parameter S_2(t) = 2Re[E_x E_y^*]', 'FontSize', 18);
grid on;
set(gca, 'FontSize', 14);
exportgraphics(gcf, 'S2_L_30um.png', 'Resolution', 300);

%% 7 Stokes 参数 S_3(t) = 2 * Im(E_x * conj(E_y)) --------

S3 = 2 * imag(E_x .* conj(E_y));
S3_norm = S3 / max(abs(S3));

figure(7);
set(gcf, 'Position', [100, 100, 1200, 400]);
%plot(t3 * 1e15, S3_norm, 'b', 'LineWidth', 2);
plot(t3 * 1e15, S3, 'b', 'LineWidth', 2);
xlabel('Time [fs]', 'FontSize', 16);
ylabel('S_3(t)', 'FontSize', 16);
title('Stokes parameter S_3(t) = 2Im[E_x E_y^*]', 'FontSize', 18);
grid on;
set(gca, 'FontSize', 14);
exportgraphics(gcf, 'S3_L_30um.png', 'Resolution', 300);

%% 8 所有 Stokes 参数叠加绘图 --------

figure(8);
set(gcf, 'Position', [100, 100, 1200, 400]);
hold on;

% plot(t3 * 1e15, S0_norm, 'k', 'LineWidth', 2, 'DisplayName', 'S_0');
% plot(t3 * 1e15, S1_norm, 'r', 'LineWidth', 2, 'DisplayName', 'S_1');
% plot(t3 * 1e15, S2_norm, 'g', 'LineWidth', 2, 'DisplayName', 'S_2');
% plot(t3 * 1e15, S3_norm, 'b', 'LineWidth', 2, 'DisplayName', 'S_3');
plot(t3 * 1e15, S0, 'k', 'LineWidth', 2, 'DisplayName', 'S_0');
plot(t3 * 1e15, S1, 'r', 'LineWidth', 2, 'DisplayName', 'S_1');
plot(t3 * 1e15, S2, 'g', 'LineWidth', 2, 'DisplayName', 'S_2');
plot(t3 * 1e15, S3, 'b', 'LineWidth', 2, 'DisplayName', 'S_3');

xlabel('Time [fs]', 'FontSize', 16);
ylabel('Stokes Parameters', 'FontSize', 16);
title('Four Stokes Parameters', 'FontSize', 18);
legend('Location', 'best', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 14);
hold off;

% 保存图像
exportgraphics(gcf, 'Stokes_All_L_30um.png', 'Resolution', 300);


