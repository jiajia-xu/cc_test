
%% 1 SF10
%Constants

omega_0 = 2.36e15;    % Central angular frequency [rad/s]
lambda = 0.8e-6;      % Wavelength [m]
Omega = 2e14;         % Spectral width [rad/s]
t0 = 1e-14;           % Pulse width [s]; T = 10fs
c = 3e8;              % Speed of light [m/s]

n = 1.7112;           % Refractive index for SF10
v_phi = c / n;        % Phase velocity [m/s]
k = omega_0 / v_phi;  % Wavevector [rad/m]
k1 = 5.83e-9;         % Group velocity [s/m]
k2 = 1.4338e-25;      % Group velocity Dispersion [s^2/m]

E0tilde = 1;
% c = 3e8;              % Speed of light [m/s]
% n = 1.453;            % Refractive index for fused silica
% v_phi = c / n;        % Phase velocity [m/s]
% omega_0 = 2.36e15;    % Central angular frequency [rad/s]
% lambda = 0.8e-6;      % Wavelength [m]
% L = 3.0e-5;           % Propagation distance [m]
% t0 = 1e-14;           % Pulse width [s]
% k = omega_0 / v_phi;  % Wavevector [rad/m]
% k1 = 5.15e-9;         % Group velocity [s/m]
% k2 = 1.99e-23;        % Group velocity Dispersion
% Omega = 2e14;         % Spectral width [rad/s]
% E0tilde = 1;
% -------- 图1: E(0,t)_x --------
t1 = linspace(-1e-13, 1e-13, 10000);  % Time array for plot 1
A0 = E0tilde * Omega * sqrt(pi);
E0_xt = A0 * exp(-t1.^2 / t0^2) .* exp(-1i * omega_0 * t1);
E0_xt_real_norm = 2 * real(E0_xt) / max(abs(2 * real(E0_xt)));  % Normalized

figure;
set(gcf, 'Position', [100, 100, 1200, 400]);  %[x y width height] → 1:3 ratio
plot(t1 * 1e15, real(E0_xt), 'b', 'LineWidth', 2);
xlabel('Time [fs]', 'FontSize', 16);
ylabel('Re[E(0,t)_x]', 'FontSize', 16);
title('Electric Field E(0,t)_x', 'FontSize', 18);
grid on;
set(gca, "FontSize", 14);
exportgraphics(gcf, 'Ex_SF10.png', 'Resolution', 300);

%% 2 no 多个L
% % -------- 图2: E_y(L,t) for multiple L values --------
% t2 = linspace(-1e-13, 4e-13, 10000);  % Time array
% L_values = [3.0e-5, 2.0e-5, 1.0e-5];  % Different propagation distances
% colors = ['r', 'g', 'b'];             % Color for each curve
% labels = ["L = 30 μm", "L = 20 μm", "L = 10 μm"];
% 
% figure;
% set(gcf, 'Position', [100, 100, 1200, 400]);  % [x y width height] → 1:3 ratio
% hold on;
% % ...（你的绘图代码保持不变）...
% for i = 1:length(L_values)
%     L = L_values(i);
%     A = 1/Omega^2 - 1i * k2 * L / 2;
%     B = k1 * L - t2;
%     Ey_Lt = E0tilde .* exp(1i * k * L) .* exp(-1i * omega_0 * t2) .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4 * A));
%     Ey_Lt_real_norm = real(Ey_Lt) / max(abs(real(Ey_Lt)));
% 
%     plot(t2 * 1e15, real(Ey_Lt), 'LineWidth', 2, 'DisplayName', labels(i)); %'Color', colors(i), 
% 
% end
% 
% xlabel('Time [fs]', 'FontSize', 16);
% ylabel('Re[E_y(L,t)]', 'FontSize', 16);
% title('Output Electric Field E_y(L,t) for Different L', 'FontSize', 18);
% legend('show', 'FontSize', 14);
% grid on;
% set(gca, 'FontSize', 14);
% hold off;
% exportgraphics(gcf, 'Ey_Lt_L_variation_SP10.png', 'Resolution', 300);


t2 = linspace(-1e-13, 40e-13, 10000);  % Time array
L_values = [10e-6, 50e-6, 200e-6, 500e-6];  % Propagation lengths
labels = ["L = 10 μm", "L = 50 μm", "L = 200 μm", "L = 500 μm"];

figure;
set(gcf, 'Position', [100, 100, 1200, 400]);
hold on;

for i = 1:length(L_values)
    L = L_values(i);
    A = 1/Omega^2 - 1i * k2 * L / 2;
    B = k1 * L - t2;
    Ey_Lt = E0tilde .* exp(1i * k * L) .* exp(-1i * omega_0 * t2) .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4 * A));
    
    %Ey_Lt_real_norm = real(Ey_Lt) / max(abs(real(Ey_Lt)));  % Normalize
    
    plot(t2 * 1e15, Ey_Lt, 'LineWidth', 2, 'DisplayName', labels(i));
end

xlabel('Time [fs]', 'FontSize', 16);
ylabel('Normalized Re[E_y(L,t)]', 'FontSize', 16);
title('E_y(L,t) in SF10 for Different Thicknesses', 'FontSize', 18);
legend('show', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 14);
hold off;
exportgraphics(gcf, 'Ey_Lt_Dispersion_SF10.png', 'Resolution', 300);


%% 2 单个L l=0.6mm
% -------- 图：E_y(L=500um, t) 显示 dispersion --------
t2 = linspace(34e-13, 36e-13, 10000);  % Time array
L = 600e-6;  % Thickness = 500 μm
label = "L = 600 μm";

figure;
set(gcf, 'Position', [100, 100, 1200, 400]);
hold on;

A = 1/Omega^2 - 1i * k2 * L / 2;
B = k1 * L - t2;
Ey_Lt = E0tilde .* exp(1i * k * L) .* exp(-1i * omega_0 * t2) ...
        .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4 * A));

%Ey_Lt_real_norm = real(Ey_Lt) / max(abs(real(Ey_Lt)));  % Normalize

plot(t2 * 1e15, Ey_Lt, 'LineWidth', 2, 'DisplayName', label);

xlabel('Time [fs]', 'FontSize', 16);
ylabel('Re[E_y(500μm,t)]', 'FontSize', 16);
title('E_y(L,t) in SF10 at L = 600 μm', 'FontSize', 18);
legend('show', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 14);
hold off;

exportgraphics(gcf, 'Ey_Lt_SF10_600um.png', 'Resolution', 300);


%% 2 单个L L=1mm
% -------- 图：E_y(L=500um, t) 显示 dispersion --------
t2 = linspace(57.3e-13, 59.3e-13, 10000);  % Time array
L = 1000e-6;  % Thickness = 500 μm
label = "L = 1 mm";

figure;
set(gcf, 'Position', [100, 100, 1200, 400]);
hold on;

A = 1/Omega^2 - 1i * k2 * L / 2;
B = k1 * L - t2;
Ey_Lt = E0tilde .* exp(1i * k * L) .* exp(-1i * omega_0 * t2) ...
        .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4 * A));

%Ey_Lt_real_norm = real(Ey_Lt) / max(abs(real(Ey_Lt)));  % Normalize

plot(t2 * 1e15, Ey_Lt, 'LineWidth', 2, 'DisplayName', label);

xlabel('Time [fs]', 'FontSize', 16);
ylabel('Re[E_y(500μm,t)]', 'FontSize', 16);
title('E_y(L,t) in SF10 at L = 1 mm', 'FontSize', 18);
legend('show', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 14);
xlim([57.3e-13* 1e15, 59.3e-13*1e15]);
hold off;

exportgraphics(gcf, 'Ey_Lt_SF10_1mm.png', 'Resolution', 300);
%% 2 单个L L=1.5mm

t2 = linspace(80.65e-13, 82.65e-13, 10000);  % Time array
L = 1400e-6;  % Thickness = 500 μm
label = "L = 1.4 mm";

figure;
set(gcf, 'Position', [100, 100, 1200, 400]);
hold on;

A = 1/Omega^2 - 1i * k2 * L / 2;
B = k1 * L - t2;
Ey_Lt = E0tilde .* exp(1i * k * L) .* exp(-1i * omega_0 * t2) ...
        .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4 * A));

%Ey_Lt_real_norm = real(Ey_Lt) / max(abs(real(Ey_Lt)));  % Normalize

plot(t2 * 1e15, Ey_Lt, 'LineWidth', 2, 'DisplayName', label);

xlabel('Time [fs]', 'FontSize', 16);
ylabel('Re[E_y(500μm,t)]', 'FontSize', 16);
title('E_y(L,t) in SF10 at L = 1.4 mm', 'FontSize', 18);
legend('show', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 14);
xlim([80.65e2, 82.65e2]);
hold off;

exportgraphics(gcf, 'Ey_Lt_SF10_1_4mm.png', 'Resolution', 300);















%% 3 - Plot S_0(t) = |E_x|^2 + |E_y|^2 for given L
% 使用 L = 1.0e-5，对应前面图一中输入 E_x 和图二中输出 E_y
% Constants
omega_0 = 2.36e15;    % Central angular frequency [rad/s]
lambda = 0.8e-6;      % Wavelength [m]
Omega = 2e14;         % Spectral width [rad/s]
t0 = 1e-14;           % Pulse width [s]; T = 10fs
c = 3e8;              % Speed of light [m/s]
L = 1.0e-3;           % Propagation distance [m]

n = 1.7112;           % Refractive index for SF10
v_phi = c / n;        % Phase velocity [m/s]
k = omega_0 / v_phi;  % Wavevector [rad/m]
k1 = 5.83e-9;         % Group velocity [s/m]
k2 = 1.4338e-25;      % Group velocity Dispersion [s^2/m]

E0tilde = 1;

A0 = E0tilde * Omega * sqrt(pi);  %below Eq.(11)

% 时间轴：我们统一选最长时间范围以便兼容
t3 = linspace(-1e-13, 7e-12, 10000);

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
figure;
set(gcf, 'Position', [100, 100, 1200, 400]);  % 1:3 横向图像比例
%plot(t3 * 1e15, S0_norm, 'm', 'LineWidth', 2);
plot(t3 * 1e15, S0, 'm', 'LineWidth', 2);
xlabel('Time [fs]', 'FontSize', 16);
ylabel('S_0(t)', 'FontSize', 16);
title('Stokes parameter S_0(t) = |E_x|^2 + |E_y|^2', 'FontSize', 18);
grid on;
set(gca, 'FontSize', 14);

% 保存图像
exportgraphics(gcf, 'S0_L_1mm.png', 'Resolution', 300);


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
exportgraphics(gcf, 'S1_L_1mm.png', 'Resolution', 300);

%% 6  Stokes 参数 S_2(t) = 2 * Re(E_x * conj(E_y)) --------
%S2 = real(conj(E_x) .* E_y) + real(conj(E_y) .* E_x);
S2 = 2 * real(E_x .* conj(E_y));
S2_norm = S2 / max(abs(S2));

figure(61);
set(gcf, 'Position', [100, 100, 1200, 400]);
%plot(t3 * 1e15, S2_norm, 'g', 'LineWidth', 2);
plot(t3 * 1e15, S2, 'g', 'LineWidth', 2);
xlabel('Time [fs]', 'FontSize', 16);
ylabel('S_2(t)', 'FontSize', 16);
title('Stokes parameter S_2(t) = 2Re[E_x E_y^*]', 'FontSize', 18);
grid on;
set(gca, 'FontSize', 14);
exportgraphics(gcf, 'S2_L_1mm.png', 'Resolution', 300);

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
exportgraphics(gcf, 'S3_L_1mm.png', 'Resolution', 300);

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
exportgraphics(gcf, 'Stokes_All_L_1mm.png', 'Resolution', 300);

%%



