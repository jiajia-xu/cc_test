%% 4 - Plot S_0(t) = |E_x|^2 + |E_y|^2 for given L
% 使用 L = 1.0e-5，对应前面图一中输入 E_x 和图二中输出 E_y
% Constants
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.5108;
L = 0.4e-5;


% Dispersion parameters
dn_dlambda = -0.0198418e6;        % /m
d2n_dlambda2 = 0.0492482e12;      % /m^2
v_phi = c / n;        % Phase velocity [m/s]
k = omega0 / v_phi;  % Wavevector [rad/m] 
k1 = (1/c) * (n + lambda0 * dn_dlambda); %group velocity
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2; % group velocity delay

% Pulse parameters
t0 = 2e-14;       %T = 20fs Figure1(b)
Omega = 2 / t0;   %below eq.11
E0tilde = 1;


% 时间轴：我们统一选最长时间范围以便兼容
t3 = linspace(-3*t0, 7*t0, 40000);

% E_x(t)
A0 = E0tilde * Omega * sqrt(pi); %below eq.11
Ex = A0 * exp(-t3.^2 / t0^2) .* exp(-1i * omega0 * t3) .* exp(1i * k0 * L);

% E_y(t) for L = 1.0e-5（与图2保持一致）
A = 1 / Omega^2 - 1i * k2 * L / 2;
B = k1 * L - t3;
Ey = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3) .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4*A));


% 计算 S0(t)
S0 = abs(Ex).^2 + abs(Ey).^2;

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
exportgraphics(gcf, 'S0_L_10um.png', 'Resolution', 300);


%% 5 Stokes 参数 S_1(t) = |E_x|^2 - |E_y|^2 --------

% 计算 S1(t)
S1 = abs(Ex).^2 - abs(Ey).^2;

% 归一化到最大绝对值（保持正负对称）
%S1_norm = S1 / max(abs(S1));

% 绘图
figure;
set(gcf, 'Position', [100, 100, 1200, 400]);  % 横向图像比例 1:3
plot(t3 * 1e15, S1, 'c', 'LineWidth', 2);
%plot(t3 * 1e15, S1_norm, 'c', 'LineWidth', 2);
xlabel('Time [fs]', 'FontSize', 16);
ylabel('S_1(t)', 'FontSize', 16);
title('Stokes parameter S_1(t) = |E_x|^2 - |E_y|^2', 'FontSize', 18);
grid on;
set(gca, 'FontSize', 14);

% 保存图像
exportgraphics(gcf, 'S1_L_10um.png', 'Resolution', 300);

%% 6  Stokes 参数 S_2(t) = 2 * Re(E_x * conj(E_y)) --------
%S2 = real(conj(E_x) .* E_y) + real(conj(E_y) .* E_x);
S2 = 2 * real(Ex .* conj(Ey));
S2_norm = S2 / max(abs(S2));

figure;
set(gcf, 'Position', [100, 100, 1200, 400]);
%plot(t3 * 1e15, S2_norm, 'g', 'LineWidth', 2);
plot(t3 * 1e15, S2, 'g', 'LineWidth', 2);
xlabel('Time [fs]', 'FontSize', 16);
ylabel('S_2(t)', 'FontSize', 16);
title('Stokes parameter S_2(t) = 2Re[E_x E_y^*]', 'FontSize', 18);
grid on;
set(gca, 'FontSize', 14);
exportgraphics(gcf, 'S2_L_10um.png', 'Resolution', 300);

%% 7 Stokes 参数 S_3(t) = 2 * Im(E_x * conj(E_y)) --------

S3 = 2 * imag(Ex .* conj(Ey));
S3_norm = S3 / max(abs(S3));

figure;
set(gcf, 'Position', [100, 100, 1200, 400]);
%plot(t3 * 1e15, S3_norm, 'b', 'LineWidth', 2);
plot(t3 * 1e15, S3, 'b', 'LineWidth', 2);
xlabel('Time [fs]', 'FontSize', 16);
ylabel('S_3(t)', 'FontSize', 16);
title('Stokes parameter S_3(t) = 2Im[E_x E_y^*]', 'FontSize', 18);
grid on;
set(gca, 'FontSize', 14);
exportgraphics(gcf, 'S3_L_10um.png', 'Resolution', 300);

%% 8 所有 Stokes 参数叠加绘图 --------

figure;
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
exportgraphics(gcf, 'Stokes_All_L_10um.png', 'Resolution', 300);


%% 9 WRONG - Visualize normalized Stokes vectors on the Poincaré sphere for L ∈ [0, 1 mm]


% E_x is fixed (independent of L)
%A0 = E0tilde * Omega * sqrt(pi); %below eq.11
Ex = A0 * exp(-t3.^2 / t0^2) .* exp(-1i * omega0 * t3) .* exp(1i * k0 * L);

% Range of L values
L_values = linspace(0, 1e-5, 200);  % From 0 to 10 um, 200 steps

% Initialize array to store normalized Stokes vectors
Stokes_normalized = zeros(length(L_values), 3);  % [S1/S0, S2/S0, S3/S0]

for idx = 1:length(L_values)
    L = L_values(idx);
    
    % E_y(t) at this L
    A = 1/Omega^2 - 1i * k2 * L / 2;
    B = k1 * L - t3;
    Ey = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3) ...
        .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4 * A));
    
    % Stokes parameters
    S0 = abs(Ex).^2 + abs(Ey).^2;
    S1 = abs(Ex).^2 - abs(Ey).^2;
    S2 = 2 * real(Ex .* conj(Ey));
    S3 = 2 * imag(Ex .* conj(Ey));
    
    % Time-averaged normalized Stokes vector
    Stokes_normalized(idx, 1) = mean(S1) / mean(S0);
    Stokes_normalized(idx, 2) = mean(S2) / mean(S0);
    Stokes_normalized(idx, 3) = mean(S3) / mean(S0);
   % Stokes_normalized(idx, 4) = mean(S0) / mean(S0);
end

% 添加标题行
header = {'S1_over_S0', 'S2_over_S0', 'S3_over_S0'};%, 'S0_over_S0'
% 保存为 table 并写入 CSV
T = array2table(Stokes_normalized, 'VariableNames', header);
writetable(T, 'stokes_normalized_1e-5.csv');

% Plotting the Poincaré sphere
figure;
[xs, ys, zs] = sphere(100);
surf(xs, ys, zs, 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % unit sphere
colormap gray;
hold on;

% Plot normalized Stokes vectors
plot3(Stokes_normalized(:,1), Stokes_normalized(:,2), Stokes_normalized(:,3), ...
    'r', 'LineWidth', 2);
scatter3(Stokes_normalized(:,1), Stokes_normalized(:,2), Stokes_normalized(:,3), ...
    20, linspace(0,1,length(L_values)), 'filled');

% Axes and labels
axis equal;
xlabel('S_1 / S_0');
ylabel('S_2 / S_0');
zlabel('S_3 / S_0');
title('Trajectory of Polarization State on the Poincaré Sphere (L=0~10um)', 'FontSize', 16);
grid on;
view(135, 30);
colorbar;
caxis([0, 1]);  % Maps to L ∈ [0, 1 mm]

% Save figure
exportgraphics(gcf, 'PoincareSphere_Trajectory.png', 'Resolution', 300);
%% 10 PoincareSphere_Trajectory
% Constants
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.5108;
L = 15e-6;

% Dispersion parameters
dn_dlambda = -0.0198e6;        % /m
d2n_dlambda2 = 0.0491e12;      % /m^2
v_phi = c / n;        
k = omega0 / v_phi;  
k1 = (1/c) * (n + lambda0 * dn_dlambda); 
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2; 

% Pulse parameters
t0 = 2e-14; 
Omega = 2 / t0; 
E0tilde = 1;

% Time axis
t3_values = linspace(-3*t0, 7*t0, 100);

% A0 is needed
A0 = E0tilde * Omega * sqrt(pi); 

% E_x is fixed
Ex = A0 * exp(-t3_values.^2 / t0^2) .* exp(-1i * omega0 * t3_values) .* exp(1i * k0 * L);

% Initialize array to store Stokes vectors
Stokes = zeros(length(t3_values), 4);  

% Compute Stokes parameters over time
for idx = 1:length(t3_values)
    t3 = t3_values(idx);
    
    A = 1/Omega^2 - 1i * k2 * L / 2;
    B = k1 * L - t3;
    Ey = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3) ...
        .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4 * A));
    
    S0 = abs(Ex(idx)).^2 + abs(Ey).^2;
    S1 = abs(Ex(idx)).^2 - abs(Ey).^2;
    S2 = 2 * real(Ex(idx) .* conj(Ey));
    S3 = 2 * imag(Ex(idx) .* conj(Ey));

    
    Stokes(idx, :) = [S1/S0, S2/S0, S3/S0, (S1^2+S2^2+S3^2)/S0^2];  % normalized
end

% Save to CSV
header = {'S1_over_S0', 'S2_over_S0', 'S3_over_S0', 'S0_over_S0'};
T = array2table(Stokes, 'VariableNames', header);
writetable(T, 'stokes_1e-5.csv');

% Plotting the Poincaré sphere trajectory using computed points
figure;
% Optional: Draw unit sphere
[xs, ys, zs] = sphere(100);
surf(xs, ys, zs, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
colormap gray;
hold on;
% 添加穿过球心的参考轴线（x, y, z）
% 添加穿过球心的参考轴线（x, y, z）为灰色实线
plot3([-1 1], [0 0], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);  % x-axis
plot3([0 0], [-1 1], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);  % y-axis
plot3([0 0], [0 0], [-1 1], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);  % z-axis


% Plot trajectory
plot3(Stokes(:,1), Stokes(:,2), Stokes(:,3), 'r', 'LineWidth', 2);
scatter3(Stokes(:,1), Stokes(:,2), Stokes(:,3), ...
    20, linspace(0,1,length(t3_values)), 'filled');

% Axes and labels
axis equal; 
xlabel('$s_1$', 'FontSize', 24, 'Interpreter', 'latex');
ylabel('$s_2$', 'FontSize', 24, 'Interpreter', 'latex');
zlabel('$s_3$', 'FontSize', 24, 'Interpreter', 'latex');
title('Trajectory of Polarization State on the Poincaré Sphere (L=15um)', 'FontSize', 16);
grid on;
view(135, 20);
%colorbar;
caxis([0, 1]);  % Maps time progression
%colorbar.Label.String = 'Normalized Time';

% Save figure
%exportgraphics(gcf, 'PoincareSphere_Trajectory_15um.png', 'Resolution', 300);
%% 10 colors 8
% Constants
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.5108;
L = 15e-6;

% Dispersion parameters
dn_dlambda = -0.0198e6;        % /m
d2n_dlambda2 = 0.0491e12;      % /m^2
v_phi = c / n;        
k = omega0 / v_phi;  
k1 = (1/c) * (n + lambda0 * dn_dlambda); 
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2; 

% Pulse parameters
t0 = 2e-14; 
Omega = 2 / t0; 
E0tilde = 1;

% Time axis
t3_values = linspace(-3*t0, 7*t0, 100);

% A0 is needed
A0 = E0tilde * Omega * sqrt(pi); 

% E_x is fixed
Ex = A0 * exp(-t3_values.^2 / t0^2) .* exp(-1i * omega0 * t3_values) .* exp(1i * k0 * L);

% Initialize array to store Stokes vectors
Stokes = zeros(length(t3_values), 4);  

% Compute Stokes parameters over time
for idx = 1:length(t3_values)
    t3 = t3_values(idx);
    
    A = 1/Omega^2 - 1i * k2 * L / 2;
    B = k1 * L - t3;
    Ey = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3) ...
        .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4 * A));
    
    S0 = abs(Ex(idx)).^2 + abs(Ey).^2;
    S1 = abs(Ex(idx)).^2 - abs(Ey).^2;
    S2 = 2 * real(Ex(idx) .* conj(Ey));
    S3 = 2 * imag(Ex(idx) .* conj(Ey));

    
    Stokes(idx, :) = [S1/S0, S2/S0, S3/S0, (S1^2+S2^2+S3^2)/S0^2];  % normalized
end

% Save to CSV
header = {'S1_over_S0', 'S2_over_S0', 'S3_over_S0', 'S0_over_S0'};
T = array2table(Stokes, 'VariableNames', header);
writetable(T, 'stokes_1e-5.csv');

% Plotting the Poincaré sphere trajectory using computed points
figure;
% Optional: Draw unit sphere
% 创建单位球
[xs, ys, zs] = sphere(100);

% 定义颜色映射（8种 RGB 颜色）
colorMap = [
    1 0 0;    % Octant 1: + + +
    0 1 0;    % Octant 2: - + +
    0 0 1;    % Octant 3: - - +
    1 1 0;    % Octant 4: + - +
    1 0 1;    % Octant 5: + + -
    0 1 1;    % Octant 6: - + -
    0.5 0.5 0.5;  % Octant 7: - - -
    1 0.5 0;      % Octant 8: + - -
];

% 初始化颜色索引矩阵
C = zeros(size(xs,1), size(xs,2), 3);

for i = 1:size(xs,1)
    for j = 1:size(xs,2)
        x = xs(i,j);
        y = ys(i,j);
        z = zs(i,j);
        
        % 确定象限编号（1~8）
        idx = 1 + (x<0) + 2*(y<0) + 4*(z<0);  % 每个位表示正负
        
        % 设置颜色
        C(i,j,:) = colorMap(idx, :);
    end
end

% 绘制彩色分区球面
surf(xs, ys, zs, C, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
hold on;

% 添加穿过球心的参考轴线（x, y, z）
% 添加穿过球心的参考轴线（x, y, z）为灰色实线
plot3([-1 1], [0 0], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);  % x-axis
plot3([0 0], [-1 1], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);  % y-axis
plot3([0 0], [0 0], [-1 1], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);  % z-axis


% Plot trajectory
plot3(Stokes(:,1), Stokes(:,2), Stokes(:,3), 'r', 'LineWidth', 2);
scatter3(Stokes(:,1), Stokes(:,2), Stokes(:,3), ...
    20, linspace(0,1,length(t3_values)), 'filled');

% Axes and labels
axis equal;
xlabel('S_1 / S_0');
ylabel('S_2 / S_0');
zlabel('S_3 / S_0');
title('Trajectory of Polarization State on the Poincaré Sphere (L=15um)', 'FontSize', 16);
grid on;
view(135, 20);
%colorbar;
caxis([0, 1]);  % Maps time progression
%colorbar.Label.String = 'Normalized Time';

% Save figure
exportgraphics(gcf, 'PoincareSphere_Trajectory_15um.png', 'Resolution', 300);

%% 10 plane ， PoincareSphere_Trajectory
% Constants
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.5108;
L = 15e-6;

% Dispersion parameters
dn_dlambda = -0.0198e6;        % /m
d2n_dlambda2 = 0.0491e12;      % /m^2
v_phi = c / n;        
k = omega0 / v_phi;  
k1 = (1/c) * (n + lambda0 * dn_dlambda); 
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2; 

% Pulse parameters
t0 = 2e-14; 
Omega = 2 / t0; 
E0tilde = 1;

% Time axis
t3_values = linspace(-3*t0, 7*t0, 100);

% A0 is needed
A0 = E0tilde * Omega * sqrt(pi); 

% E_x is fixed
Ex = A0 * exp(-t3_values.^2 / t0^2) .* exp(-1i * omega0 * t3_values) .* exp(1i * k0 * L);

% Initialize array to store Stokes vectors
Stokes = zeros(length(t3_values), 4);  

% Compute Stokes parameters over time
for idx = 1:length(t3_values)
    t3 = t3_values(idx);
    
    A = 1/Omega^2 - 1i * k2 * L / 2;
    B = k1 * L - t3;
    Ey = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3) ...
        .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4 * A));
    
    S0 = abs(Ex(idx)).^2 + abs(Ey).^2;
    S1 = abs(Ex(idx)).^2 - abs(Ey).^2;
    S2 = 2 * real(Ex(idx) .* conj(Ey));
    S3 = - 2 * imag(Ex(idx) .* conj(Ey));

    
    Stokes(idx, :) = [S1/S0, S2/S0, S3/S0, (S1^2+S2^2+S3^2)/S0^2];  % normalized
end

% Save to CSV
header = {'S1_over_S0', 'S2_over_S0', 'S3_over_S0', 'S0_over_S0'};
T = array2table(Stokes, 'VariableNames', header);
writetable(T, 'stokes_1e-5.csv');

% Plotting the Poincaré sphere trajectory using computed points
figure;
% Optional: Draw unit sphere
[xs, ys, zs] = sphere(100);
surf(xs, ys, zs, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
colormap gray;
hold on;

%
% 在单位球内画三个浅蓝色的参考面 (x=0, y=0, z=0)

% Z=0 平面 (XY 平面)
patch([-1 1 1 -1], [-1 -1 1 1], [0 0 0 0], ...
    [0.7 0.85 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Y=0 平面 (XZ 平面)
patch([-1 1 1 -1], [0 0 0 0], [-1 -1 1 1], ...
    [0.7 0.85 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% X=0 平面 (YZ 平面)
patch([0 0 0 0], [-1 1 1 -1], [-1 -1 1 1], ...
    [0.7 0.85 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

%

% Plot trajectory
plot3(Stokes(:,1), Stokes(:,2), Stokes(:,3), 'r', 'LineWidth', 2);
scatter3(Stokes(:,1), Stokes(:,2), Stokes(:,3), ...
    20, linspace(0,1,length(t3_values)), 'filled');

% Axes and labels
axis equal;
xlabel('S_1 / S_0');
ylabel('S_2 / S_0');
zlabel('S_3 / S_0');
title('Trajectory of Polarization State on the Poincaré Sphere (L=15um)', 'FontSize', 16);
grid on;
view(135, 20);
%colorbar;
caxis([0, 1]);  % Maps time progression
%colorbar.Label.String = 'Normalized Time';

% Save figure
exportgraphics(gcf, 'PoincareSphere_Trajectory_15um.png', 'Resolution', 300);
%% 10Good Four colors
% --------------------- Constants ---------------------
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.5108;
L = 15e-6;

% ---------------- Dispersion Parameters --------------
dn_dlambda = -0.0198e6;        % /m
d2n_dlambda2 = 0.0491e12;      % /m^2
v_phi = c / n;        
k = omega0 / v_phi;  
k1 = (1/c) * (n + lambda0 * dn_dlambda); 
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2; 

% ------------------ Pulse Parameters ------------------
t0 = 2e-14; 
Omega = 2 / t0; 
E0tilde = 1;
A0 = E0tilde * Omega * sqrt(pi); 

% ------------------ Time Axis ------------------
t3_values = linspace(-3*t0, 7*t0, 100);
Ex = A0 * exp(-t3_values.^2 / t0^2) .* exp(-1i * omega0 * t3_values) .* exp(1i * k0 * L);

% ------------------ Stokes Vector Calculation ------------------
Stokes = zeros(length(t3_values), 4);  

for idx = 1:length(t3_values)
    t3 = t3_values(idx);
    
    A = 1/Omega^2 - 1i * k2 * L / 2;
    B = k1 * L - t3;
    Ey = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3) ...
        .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4 * A));
    
    S0 = abs(Ex(idx)).^2 + abs(Ey).^2;
    S1 = abs(Ex(idx)).^2 - abs(Ey).^2;
    S2 = 2 * real(Ex(idx) .* conj(Ey));
    S3 = -2 * imag(Ex(idx) .* conj(Ey));
    
    Stokes(idx, :) = [S1/S0, S2/S0, S3/S0, (S1^2+S2^2+S3^2)/S0^2];  % normalized
end

% ------------------ Save CSV ------------------
header = {'S1_over_S0', 'S2_over_S0', 'S3_over_S0', 'S0_over_S0'};
T = array2table(Stokes, 'VariableNames', header);
writetable(T, 'stokes_1e-5.csv');

% ------------------ Plotting ------------------
figure;

% ---- Create colored sphere surface (4 quadrants) ----
[xs, ys, zs] = sphere(100);

% Define 4 colors for 4 x-y quadrants
colors4 = [
    1 0.6 0.6;  % Q1: x>0, y>0
    0.6 0.6 1;  % Q2: x<0, y>0
    0.6 1 0.6;  % Q3: x<0, y<0
    1 1 0.6     % Q4: x>0, y<0
];

% Construct RGB color matrix
C = zeros([size(xs), 3]);

for i = 1:size(xs,1)
    for j = 1:size(xs,2)
        x = xs(i,j);
        y = ys(i,j);

        if x >= 0 && y >= 0
            C(i,j,:) = colors4(1,:);
        elseif x < 0 && y >= 0
            C(i,j,:) = colors4(2,:);
        elseif x < 0 && y < 0
            C(i,j,:) = colors4(3,:);
        else
            C(i,j,:) = colors4(4,:);
        end
    end
end

% Plot colored unit sphere
surf(xs, ys, zs, C, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
hold on;

% ---- Add center axis lines ----
plot3([-1 1], [0 0], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);  % x
plot3([0 0], [-1 1], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);  % y
plot3([0 0], [0 0], [-1 1], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);  % z

% ---- Plot Stokes trajectory ----
plot3(Stokes(:,1), Stokes(:,2), Stokes(:,3), 'r', 'LineWidth', 2);
scatter3(Stokes(:,1), Stokes(:,2), Stokes(:,3), ...
        20, linspace(0,1,length(t3_values)), 'filled');

% ---- Axes and view settings ----
axis equal;
set(gca, 'FontSize', 14);  %'FontName', 'Times New Roman',  
xlabel('$s_1$', 'FontSize', 24, 'Interpreter', 'latex');
ylabel('$s_2$', 'FontSize', 24, 'Interpreter', 'latex');
zlabel('$s_3$', 'FontSize', 24, 'Interpreter', 'latex');
%title('Trajectory of Polarization State on the Poincaré Sphere (L = 5 μm)', 'FontSize', 14);
grid on;
view(135, 60);
colorbar;
caxis([0, 1]);
colorbar.Label.String = 'Normalized Time';

% ---- Save figure ----
%exportgraphics(gcf, 'PoincareSphere_Trajectory_color15um.png', 'Resolution', 300);

%% 10Good 3L 4colors trajectory together 
% --------------------- Constants ---------------------
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.51078;

% ---------------- Dispersion Parameters --------------
dn_dlambda = -0.0198418e6;        % /m
d2n_dlambda2 = 0.0492482e12;      % /m^2
v_phi = c / n;
k = omega0 / v_phi;  
k1 = (1/c) * (n + lambda0 * dn_dlambda); 
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2; 

% ------------------ Pulse Parameters ------------------
t0 = 2e-14; 
Omega = 2 / t0; 
%E0tilde = 1;
%A0 = E0tilde * Omega * sqrt(pi); 
A0 = 1/2;
E0tilde = A0/(Omega * sqrt(pi));
t3_values = linspace(-3*t0, 7*t0, 201);

% ------------------ Plot Setup ------------------
figure;

% ---- Create colored unit sphere surface ----
[xs, ys, zs] = sphere(100);
colors4 = [1 0.6 0.6; 0.6 0.6 1; 0.6 1 0.6; 1 1 0.6];
C = zeros([size(ys), 3]);
for i = 1:size(ys,1)
    for j = 1:size(ys,2)
         y = ys(i,j); z = zs(i,j);
        if y >= 0 && z >= 0
            C(i,j,:) = colors4(1,:);
        elseif y < 0 && z >= 0
            C(i,j,:) = colors4(2,:);
        elseif y < 0 && z < 0
            C(i,j,:) = colors4(3,:);
        else
            C(i,j,:) = colors4(4,:);
        end
    end
end
ax = axes; hold(ax,'on');

% ---- Create colored unit sphere surface ----
surf(xs, ys, zs, C, 'EdgeColor','none', 'FaceAlpha',0.15, ...
     'HandleVisibility','off');  % 不进入 legend

% ---- Add axis lines ----
plot3([-1 1],[0 0],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, ...
      'HandleVisibility','off');
plot3([0 0],[-1 1],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, ...
      'HandleVisibility','off');
plot3([0 0],[0 0],[-1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, ...
      'HandleVisibility','off');

% ---- Loop over different L values ----
L_list = [5e-6, 10e-6, 15e-6];
col = [0 1 1;   % cyan
       0 0 1;   % blue
       1 0 0];  % red

hLines = gobjects(numel(L_list),1);   % 保存三条线的句柄
legend_entries = cell(1,numel(L_list));

for li = 1:numel(L_list)
    L = L_list(li);
    Ex = A0 * exp(-t3_values.^2/t0^2) .* exp(-1i*omega0*t3_values) .* exp(1i*k0*L);
    Stokes = zeros(numel(t3_values),4);

    for idx = 1:numel(t3_values)
        t3 = t3_values(idx);
        A = 1/Omega^2 - 1i * k2 * L / 2;
        B = k1 * L - t3;
        Ey = E0tilde .* exp(1i*k*L) .* exp(-1i*omega0*t3) ...
            .* sqrt(pi./A) .* exp(-B.^2./(4*A));
        S0 = abs(Ex(idx)).^2 + abs(Ey).^2;
        S1 = abs(Ex(idx)).^2 - abs(Ey).^2;
        S2 = 2*real(Ex(idx).*conj(Ey));
        S3 = -2*imag(Ex(idx).*conj(Ey));
        Stokes(idx,:) = [S1/S0, S2/S0, S3/S0, (S1^2+S2^2+S3^2)/S0^2];
    end

    % 轨迹实线（进入 legend）
    hLines(li) = plot3(Stokes(:,1),Stokes(:,2),Stokes(:,3), ...
                       'Color',col(li,:), 'LineWidth',2);

    % 彩色散点（不进入 legend）
    scatter3(Stokes(:,1),Stokes(:,2),Stokes(:,3),20, ...
             linspace(0,1,numel(t3_values)), 'filled', ...
             'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.6, ...
             'HandleVisibility','off');

    legend_entries{li} = sprintf('$L =$ %.0f $\\mu\\mathrm{m}$', L*1e6);
end
%view(120, 45);
view(120, 30);
%view(0, 90);%from top view
legend(hLines, legend_entries, 'Location','northeast', 'Box', 'on', 'Interpreter', 'latex');

% ---- Final Plot Settings ----
axis equal;
set(gca, 'FontSize', 14);  %'FontName', 'Times New Roman',  
xticks([-1, 0, 1]);
yticks([-1, 0, 1]);
zticks([-1, 0, 1]);
xlabel('$s_1$', 'FontSize', 24, 'Interpreter', 'latex');
ylabel('$s_2$', 'FontSize', 24, 'Interpreter', 'latex');
zlabel('$s_3$', 'FontSize', 24, 'Interpreter', 'latex');
%title('Polarization Trajectories on Poincaré Sphere for Different L');
grid on;

text(1.1, 0, -0.1, 'A', 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
text(-1.1, 0, 0, 'B', 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');


%exportgraphics(gcf, 'PoincareSphere_3Trajectory_color.png', 'Resolution',300)
%exportgraphics(gcf, 'PoincareSphere_3Trajectory_topview_color.png', 'Resolution',300)

%% 10GGood 3L 4colors trajectory together - spline -4 L values
% --------------------- Constants ---------------------
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.51078;

% ---------------- Dispersion Parameters --------------
dn_dlambda = -0.0198418e6;        % /m
d2n_dlambda2 = 0.0492482e12;      % /m^2
v_phi = c / n;
k = omega0 / v_phi;  
k1 = (1/c) * (n - lambda0 * dn_dlambda); 
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2; 

% ------------------ Pulse Parameters ------------------
t0 = 2e-14; 
Omega = 2 / t0; 
%E0tilde = 1;
%A0 = E0tilde * Omega * sqrt(pi); 
A0 = 1/2;
E0tilde = A0/(Omega * sqrt(pi));
t3_values = linspace(-1*t0, 3*t0, 41);

% ------------------ Plot Setup ------------------
figure;

% ---- Create colored unit sphere surface ----
[xs, ys, zs] = sphere(100);
colors4 = [1 0.6 0.6; 0.6 0.6 1; 0.6 1 0.6; 1 1 0.6];
C = zeros([size(ys), 3]);
for i = 1:size(ys,1)
    for j = 1:size(ys,2)
         y = ys(i,j); z = zs(i,j);
        if y >= 0 && z >= 0
            C(i,j,:) = colors4(1,:);
        elseif y < 0 && z >= 0
            C(i,j,:) = colors4(2,:);
        elseif y < 0 && z < 0
            C(i,j,:) = colors4(3,:);
        else
            C(i,j,:) = colors4(4,:);
        end
    end
end
ax = axes; hold(ax,'on');

% ---- Create colored unit sphere surface ----
surf(xs, ys, zs, C, 'EdgeColor','none', 'FaceAlpha',0.15, ...
     'HandleVisibility','off');  % 不进入 legend

% ---- Add axis lines ----
plot3([-1 1],[0 0],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, ...
      'HandleVisibility','off');
plot3([0 0],[-1 1],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, ...
      'HandleVisibility','off');
plot3([0 0],[0 0],[-1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, ...
      'HandleVisibility','off');

% ---- Loop over different L values ----
L_list = [5e-6, 10e-6, 15e-6, 17e-6];
col = [0 1 1;   % cyan
       0 0 1;   % blue
       1 0 0;   % red
       0.6 0 0.6];% purple 

hLines = gobjects(numel(L_list),1);     % 保存三条线的句柄
legend_entries = cell(1,numel(L_list));

for li = 1:numel(L_list)
    L = L_list(li);
    Ex = A0 * exp(-t3_values.^2/t0^2) .* exp(-1i*omega0*t3_values) .* exp(1i*k0*L);
    Stokes = zeros(numel(t3_values),4);

    for idx = 1:numel(t3_values)
        t3 = t3_values(idx);
        A = 1/Omega^2 - 1i * k2 * L / 2;
        B = k1 * L - t3;
        Bt(idx) = k1 * L - t3;
        Ey = E0tilde .* exp(1i*k*L) .* exp(-1i*omega0*t3) ...
            .* sqrt(pi./A) .* exp(-B.^2./(4*A));
        S0 = abs(Ex(idx)).^2 + abs(Ey).^2;
        S1 = abs(Ex(idx)).^2 - abs(Ey).^2;
        S2 = 2*real(Ex(idx).*conj(Ey));
        S3 = - 2*imag(Ex(idx).*conj(Ey));
        Stokes(idx,:) = [S1/S0, S2/S0, S3/S0, (S1^2+S2^2+S3^2)/S0^2];
    end

    % 使用样条插值进行平滑
    t_idx = 1:numel(t3_values);                     % 原始点参数
    tt = linspace(1, numel(t3_values), 20*numel(t3_values));  % 插值参数

    % 对每个分量进行样条插值
    s1_spline = spline(t_idx, Stokes(:,1), tt);
    s2_spline = spline(t_idx, Stokes(:,2), tt);
    s3_spline = spline(t_idx, Stokes(:,3), tt);

    % 平滑轨迹实线（进入 legend）
    hLines(li) = plot3(s1_spline, s2_spline, s3_spline, ...
                   'Color', col(li,:), 'LineWidth', 2);


    % 彩色散点（不进入 legend）
    scatter3(Stokes(:,1),Stokes(:,2),Stokes(:,3), 30, ...
             linspace(0,1,numel(t3_values)), 'filled', ...
             'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.9, ...
             'HandleVisibility','off');


    % % 彩色散点（颜色与线条匹配，不进入 legend）
    %      scatter3(Stokes(:,1), Stokes(:,2), Stokes(:,3), 20, ...
    %      repmat(col(li,:), size(Stokes,1), 1), 'filled', ...
    %      'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.6, ...
    %      'HandleVisibility','off');


    legend_entries{li} = sprintf('$L =$ %.0f $\\mu\\mathrm{m}$', L*1e6);
end
%view(120, 45);
view(120, 30);
%view(0, 90);%from top view
legend(hLines, legend_entries, 'Location','northeast', 'Box', 'on', 'Interpreter', 'latex');

% ---- Final Plot Settings ----
axis equal;
set(gca, 'FontSize', 14);  %'FontName', 'Times New Roman',  
xticks([-1, 0, 1]);
yticks([-1, 0, 1]);
zticks([-1, 0, 1]);
xlabel('$s_1$', 'FontSize', 24, 'Interpreter', 'latex');
ylabel('$s_2$', 'FontSize', 24, 'Interpreter', 'latex');
zlabel('$s_3$', 'FontSize', 24, 'Interpreter', 'latex');
%title('Polarization Trajectories on Poincaré Sphere for Different L');
grid on;

text(1.1, 0, -0.1, 'I', 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
text(-1.1, 0.1, 0, 'F', 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');


exportgraphics(gcf, 'PoincareSphere_4Trajectory_color_spline.png', 'Resolution',300)

%% 10 （251021）GGood colorbar -- 3L 4colors trajectory together - spline -4 L values
% --------------------- Constants ---------------------
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.51078;

% ---------------- Dispersion Parameters --------------
dn_dlambda = -0.0198418e6;        % /m
d2n_dlambda2 = 0.0492482e12;      % /m^2
v_phi = c / n;
k = omega0 / v_phi;  
k1 = (1/c) * (n - lambda0 * dn_dlambda); 
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2; 

% ------------------ Pulse Parameters ------------------
t0 = 2e-14; 
Omega = 2 / t0; 
A0 = 1/2;
E0tilde = A0/(Omega * sqrt(pi));

% 时间轴与用于着色的 fs 版本
t3_values = linspace(-1*t0, 3*t0, 41);   % 单位: s
t_fs = t3_values * 1e15;                 % 单位: fs，用于颜色映射

% ------------------ Plot Setup ------------------
figure;

% ---- Create colored unit sphere surface ----
[xs, ys, zs] = sphere(100);
colors4 = [1 0.6 0.6; 0.6 0.6 1; 0.6 1 0.6; 1 1 0.6];
C = zeros([size(ys), 3]);
for i = 1:size(ys,1)
    for j = 1:size(ys,2)
        y = ys(i,j); z = zs(i,j);
        if y >= 0 && z >= 0
            C(i,j,:) = colors4(1,:);
        elseif y < 0 && z >= 0
            C(i,j,:) = colors4(2,:);
        elseif y < 0 && z < 0
            C(i,j,:) = colors4(3,:);
        else
            C(i,j,:) = colors4(4,:);
        end
    end
end
ax = axes; hold(ax,'on');

% ---- Create colored unit sphere surface ----
surf(xs, ys, zs, C, 'EdgeColor','none', 'FaceAlpha',0.15, ...
     'HandleVisibility','off');  % 不进入 legend

% ---- Add axis lines ----
plot3([-1 1],[0 0],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, ...
      'HandleVisibility','off');
plot3([0 0],[-1 1],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, ...
      'HandleVisibility','off');
plot3([0 0],[0 0],[-1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, ...
      'HandleVisibility','off');

% ---- Loop over different L values ----
L_list = [5e-6, 10e-6, 15e-6, 17e-6];
col = [0 1 1;         % cyan
       0 0 1;         % blue
       1 0 0;         % red
       0.6 0 0.6];    % purple 

hLines = gobjects(numel(L_list),1);   % 保存三条线的句柄
legend_entries = cell(1,numel(L_list));

for li = 1:numel(L_list)
    L = L_list(li);
    Ex = A0 * exp(-t3_values.^2/t0^2) .* exp(-1i*omega0*t3_values) .* exp(1i*k0*L);
    Stokes = zeros(numel(t3_values),4);

    for idx = 1:numel(t3_values)
        t3 = t3_values(idx);
        A = 1/Omega^2 - 1i * k2 * L / 2;
        B = k1 * L - t3;
        Bt(idx) = k1 * L - t3; %#ok<AGROW>
        Ey = E0tilde .* exp(1i*k*L) .* exp(-1i*omega0*t3) ...
            .* sqrt(pi./A) .* exp(-B.^2./(4*A));
        S0 = abs(Ex(idx)).^2 + abs(Ey).^2;
        S1 = abs(Ex(idx)).^2 - abs(Ey).^2;
        S2 = 2*real(Ex(idx).*conj(Ey));
        S3 = - 2*imag(Ex(idx).*conj(Ey));
        Stokes(idx,:) = [S1/S0, S2/S0, S3/S0, (S1^2+S2^2+S3^2)/S0^2];
    end

    % 使用样条插值进行平滑
    t_idx = 1:numel(t3_values);                             % 原始点参数
    tt = linspace(1, numel(t3_values), 20*numel(t3_values));% 插值参数

    % 对每个分量进行样条插值
    s1_spline = spline(t_idx, Stokes(:,1), tt);
    s2_spline = spline(t_idx, Stokes(:,2), tt);
    s3_spline = spline(t_idx, Stokes(:,3), tt);

    % 平滑轨迹实线（进入 legend）
    hLines(li) = plot3(s1_spline, s2_spline, s3_spline, ...
                       'Color', col(li,:), 'LineWidth', 2);

    % 彩色散点（用时间 fs 作颜色映射，不进入 legend）
    scatter3(Stokes(:,1), Stokes(:,2), Stokes(:,3), 30, ...
             t_fs, 'filled', ...
             'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.9, ...
             'HandleVisibility','off');

    legend_entries{li} = sprintf('$L =$ %.0f $\\mu\\mathrm{m}$', L*1e6);
end

view(120, 30);
legend(hLines, legend_entries, 'Location','northeast', 'Box', 'on', 'Interpreter', 'latex');

% ---- Final Plot Settings ----
axis equal;
set(gca, 'FontSize', 14);
xticks([-1, 0, 1]); yticks([-1, 0, 1]); zticks([-1, 0, 1]);
xlabel('$s_1$', 'FontSize', 24, 'Interpreter', 'latex');
ylabel('$s_2$', 'FontSize', 24, 'Interpreter', 'latex');
zlabel('$s_3$', 'FontSize', 24, 'Interpreter', 'latex');
grid on;

text(1.1, 0, -0.1, 'I', 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
text(-1.1, 0.1, 0, 'F', 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');

% ---- 添加 colorbar 表示时间（fs）----
colormap(parula);                 % 可改为 turbo/jet 等
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = 'evolution time~[fs]';
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 14;

% 颜色条范围与刻度：精确从 -20 fs 到 60 fs
caxis([-20, 60]);
cb.Ticks = -20:10:60;

exportgraphics(gcf, 'PoincareSphere_4Trajectory_color_spline_251021.png', 'Resolution',300)

%% 10GGood 251021 3L 4colors trajectory together - spline -4 L values (dot-style colorbar)
% --------------------- Constants ---------------------
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.51078;

% ---------------- Dispersion Parameters --------------
dn_dlambda = -0.0198418e6;        % /m
d2n_dlambda2 = 0.0492482e12;      % /m^2
v_phi = c / n;
k = omega0 / v_phi;  
k1 = (1/c) * (n - lambda0 * dn_dlambda); 
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2; 

% ------------------ Pulse Parameters ------------------
t0 = 2e-14; 
Omega = 2 / t0; 
A0 = 1/2;
E0tilde = A0/(Omega * sqrt(pi));

% 时间轴与用于着色的 fs 版本
t3_values = linspace(-1*t0, 3*t0, 41);   % 单位: s
t_fs = t3_values * 1e15;                 % 单位: fs，用于颜色映射

% ------------------ Plot Setup ------------------
figure;

% ---- Create colored unit sphere surface ----
[xs, ys, zs] = sphere(100);
colors4 = [1 0.6 0.6; 0.6 0.6 1; 0.6 1 0.6; 1 1 0.6];
C = zeros([size(ys), 3]);
for i = 1:size(ys,1)
    for j = 1:size(ys,2)
        y = ys(i,j); z = zs(i,j);
        if y >= 0 && z >= 0
            C(i,j,:) = colors4(1,:);
        elseif y < 0 && z >= 0
            C(i,j,:) = colors4(2,:);
        elseif y < 0 && z < 0
            C(i,j,:) = colors4(3,:);
        else
            C(i,j,:) = colors4(4,:);
        end
    end
end
ax = axes; hold(ax,'on');

% ---- Create colored unit sphere surface ----
surf(xs, ys, zs, C, 'EdgeColor','none', 'FaceAlpha',0.15, ...
     'HandleVisibility','off');

% ---- Add axis lines ----
plot3([-1 1],[0 0],[0 0],'--','Color',[0.4 0.4 0.4],'LineWidth',0.8, ...
      'HandleVisibility','off');
plot3([0 0],[-1 1],[0 0],'--','Color',[0.4 0.4 0.4],'LineWidth',0.8, ...
      'HandleVisibility','off');
plot3([0 0],[0 0],[-1 1],'--','Color',[0.4 0.4 0.4],'LineWidth',0.8, ...
      'HandleVisibility','off');

% ---- Loop over different L values ----
L_list = [5e-6, 10e-6, 15e-6, 17e-6];
col = [0 1 1;         % cyan
       0 0 1;         % blue
       1 0 0;         % red
       0.6 0 0.6];    % purple 

hLines = gobjects(numel(L_list),1);
legend_entries = cell(1,numel(L_list));

for li = 1:numel(L_list)
    L = L_list(li);
    Ex = A0 * exp(-t3_values.^2/t0^2) .* exp(-1i*omega0*t3_values) .* exp(1i*k0*L);
    Stokes = zeros(numel(t3_values),4);

    for idx = 1:numel(t3_values)
        t3 = t3_values(idx);
        A = 1/Omega^2 - 1i * k2 * L / 2;
        B = k1 * L - t3;
        Bt(idx) = k1 * L - t3; %#ok<AGROW>
        Ey = E0tilde .* exp(1i*k*L) .* exp(-1i*omega0*t3) ...
            .* sqrt(pi./A) .* exp(-B.^2./(4*A));
        S0 = abs(Ex(idx)).^2 + abs(Ey).^2;
        S1 = abs(Ex(idx)).^2 - abs(Ey).^2;
        S2 = 2*real(Ex(idx).*conj(Ey));
        S3 = - 2*imag(Ex(idx).*conj(Ey));
        Stokes(idx,:) = [S1/S0, S2/S0, S3/S0, (S1^2+S2^2+S3^2)/S0^2];
    end

    % 样条插值
    t_idx = 1:numel(t3_values);
    tt = linspace(1, numel(t3_values), 20*numel(t3_values));
    s1_spline = spline(t_idx, Stokes(:,1), tt);
    s2_spline = spline(t_idx, Stokes(:,2), tt);
    s3_spline = spline(t_idx, Stokes(:,3), tt);

    % 平滑轨迹实线（进入 legend）
    hLines(li) = plot3(s1_spline, s2_spline, s3_spline, ...
                       'Color', col(li,:), 'LineWidth', 2);

    % 彩色散点（用时间 fs 作颜色映射，不进入 legend）
    scatter3(Stokes(:,1), Stokes(:,2), Stokes(:,3), 30, ...
             t_fs, 'filled', ...
             'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.9, ...
             'HandleVisibility','off');

    legend_entries{li} = sprintf('$L =$ %.0f $\\mu\\mathrm{m}$', L*1e6);
end

% 视图 & 图例
view(120, 30);
legend(hLines, legend_entries, 'Location','northeast', 'Box', 'on', 'Interpreter', 'latex');

% ---- Final Plot Settings ----
axis equal;
set(gca, 'FontSize', 15);
xticks([-1, 0, 1]); yticks([-1, 0, 1]); zticks([-1, 0, 1]);
xlabel('$s_1$', 'FontSize', 24, 'Interpreter', 'latex');
ylabel('$s_2$', 'FontSize', 24, 'Interpreter', 'latex');
zlabel('$s_3$', 'FontSize', 24, 'Interpreter', 'latex');

ax = gca;
grid on
ax.Color = 'w';                % 背景白
ax.GridColor = [0.6 0.6 0.6];  % 深灰网格
ax.GridAlpha = 0.6;            % 半透明
ax.LineWidth = 1;              % 加粗坐标轴

text(1.1, 0, -0.1, 'I', 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
text(-1.1, 0.1, 0, 'F', 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');

% ---- 颜色映射设置：确保主图颜色与“点状 colorbar”一致 ----
cmap = parula(256);           % 也可改为 turbo(256)
colormap(cmap);
caxis([-20, 60]);             % fs 范围固定为 -20 到 60

% ---- 点状 colorbar：在右侧绘制彩色点 ----
t_ticks = [-20, -10, 0, 10, 20, 30, 40, 50, 60];        % 显示刻度不变
axpos = get(gca,'Position');           % 主轴位置

% 调整位置 & 大小（更紧凑）
dot_ax = axes('Position',[axpos(1)+axpos(3)+0.05, axpos(2)+0.04*axpos(4), 0.035, 0.8*axpos(4)]);
hold(dot_ax,'on');

% 依据相同 colormap & caxis 生成对应颜色
t_all = linspace(-20, 60, size(cmap,1));

% 🔸 更密集的 dot 时间点（例如每 2 fs 一个）
t_dense = linspace(-20, 60, 41);
dot_dense_colors = interp1(t_all, cmap, t_dense, 'linear', 'extrap');

% 绘制密集的小点（颜色平滑渐变）
scatter(dot_ax, zeros(size(t_dense)), t_dense, 36, dot_dense_colors, 'filled', 'MarkerEdgeColor','none');

% 🔸 叠加主要刻度点（略大一些，用于视觉参考）
dot_colors = interp1(t_all, cmap, t_ticks, 'linear', 'extrap');
scatter(dot_ax, zeros(size(t_ticks)), t_ticks, 50, dot_colors, 'filled', 'MarkerEdgeColor','k');

% 坐标轴设置
xlim(dot_ax, [-1 1]); ylim(dot_ax, [-25 65]);
dot_ax.YTick = t_ticks;
dot_ax.YTickLabel = arrayfun(@(x)sprintf('%d',x), t_ticks, 'UniformOutput', false);
dot_ax.YLabel.String = 'time~[fs]';
dot_ax.YLabel.Interpreter = 'latex';
dot_ax.YLabel.FontSize = 14;
dot_ax.XColor = 'none'; 
dot_ax.Box = 'off';
dot_ax.YDir = 'normal';
%dot_ax.YAxisLocation = 'right';       % 👈 新增：让标签和刻度都在右侧
set(dot_ax, 'TickLabelInterpreter','latex', 'FontSize', 14);


% （可选）标题或说明
% title(dot_ax, 'time', 'Interpreter','latex');

% ---- 导出 ----
exportgraphics(gcf, 'PoincareSphere_4Trajectory_color_spline_dotsLegend.png', 'Resolution', 300);


%% 10 Test2.8um-4.3um L-16lines

clc; clear;

% ---------------- 基本参数 ----------------
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.51078;

dn_dlambda = -0.0198418e6;        % /m
d2n_dlambda2 = 0.0492482e12;      % /m^2
v_phi = c / n;
k = omega0 / v_phi;  
k1 = (1/c) * (n - lambda0 * dn_dlambda); 
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2; 

t0 = 2e-14; 
Omega = 2 / t0; 
A0 = 1/2;
E0tilde = A0 / (Omega * sqrt(pi));
t3_values = linspace(-3*t0, 7*t0, 201);

% ---------------- Poincaré 球面设置 ----------------
figure;
[xs, ys, zs] = sphere(100);
colors4 = [1 0.6 0.6; 0.6 0.6 1; 0.6 1 0.6; 1 1 0.6];
C = zeros([size(ys), 3]);
for i = 1:size(ys,1)
    for j = 1:size(ys,2)
         y = ys(i,j); z = zs(i,j);
        if y >= 0 && z >= 0
            C(i,j,:) = colors4(1,:);
        elseif y < 0 && z >= 0
            C(i,j,:) = colors4(2,:);
        elseif y < 0 && z < 0
            C(i,j,:) = colors4(3,:);
        else
            C(i,j,:) = colors4(4,:);
        end
    end
end
ax = axes; hold(ax,'on');
surf(xs, ys, zs, C, 'EdgeColor','none', 'FaceAlpha',0.15, ...
     'HandleVisibility','off');

% ---- 轴线 ----
plot3([-1 1],[0 0],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, 'HandleVisibility','off');
plot3([0 0],[-1 1],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, 'HandleVisibility','off');
plot3([0 0],[0 0],[-1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, 'HandleVisibility','off');

% ---------------- 主循环：11 条轨迹 ----------------
L_list = linspace(2.8e-6, 4.3e-6, 16);
colors = jet(numel(L_list));
hLines = gobjects(numel(L_list),1);
legend_entries = cell(1,numel(L_list));


for li = 1:numel(L_list)
    L = L_list(li);
    Ex = A0 * exp(-t3_values.^2/t0^2) .* exp(-1i*omega0*t3_values) .* exp(1i*k0*L);
    Stokes = zeros(numel(t3_values),4);

    for idx = 1:numel(t3_values)
        t3 = t3_values(idx);
        A = 1/Omega^2 - 1i * k2 * L / 2;
        B = k1 * L - t3;
        Ey = E0tilde .* exp(1i*k*L) .* exp(-1i*omega0*t3) ...
            .* sqrt(pi./A) .* exp(-B.^2./(4*A));
        S0 = abs(Ex(idx)).^2 + abs(Ey).^2;
        S1 = abs(Ex(idx)).^2 - abs(Ey).^2;
        S2 = 2*real(Ex(idx).*conj(Ey));
        S3 = - 2*imag(Ex(idx).*conj(Ey));
        Stokes(idx,:) = [S1/S0, S2/S0, S3/S0, (S1^2+S2^2+S3^2)/S0^2];
    end

    % 插值平滑轨迹
    t_idx = 1:numel(t3_values);
    tt = linspace(1, numel(t3_values), 20*numel(t3_values));
    s1_spline = spline(t_idx, Stokes(:,1), tt);
    s2_spline = spline(t_idx, Stokes(:,2), tt);
    s3_spline = spline(t_idx, Stokes(:,3), tt);

    % 平滑轨迹线
    hLines(li) = plot3(s1_spline, s2_spline, s3_spline, ...
                       'Color', colors(li,:), 'LineWidth', 1.5);

    % 散点轨迹（不进 legend）
    scatter3(Stokes(:,1),Stokes(:,2),Stokes(:,3),20, ...
             linspace(0,1,numel(t3_values)), 'filled', ...
             'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.6, ...
             'HandleVisibility','off');

    legend_entries{li} = sprintf('$L = %.1f\\,\\mu\\mathrm{m}$', L*1e6);
    % --- 添加每条轨迹的中心点标注 ---
mid_idx = round(numel(t3_values)/3);
s1c = Stokes(mid_idx,1);
s2c = Stokes(mid_idx,2);
s3c = Stokes(mid_idx,3);
text(s1c, s2c, s3c, ...
     sprintf('%.1f\\,$\\mu$m', L*1e6), ...
     'FontSize', 10, 'Interpreter','latex', ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','middle', ...
     'Color', colors(li,:));

end

% ---------------- 图像设置 ----------------
view(120, 30);
%legend(hLines, legend_entries, 'Location','northeast', 'Box','on', 'Interpreter','latex');

axis equal;
set(gca, 'FontSize', 14);
xticks([-1, 0, 1]); yticks([-1, 0, 1]); zticks([-1, 0, 1]);
xlabel('$s_1$', 'FontSize', 24, 'Interpreter','latex');
ylabel('$s_2$', 'FontSize', 24, 'Interpreter','latex');
zlabel('$s_3$', 'FontSize', 24, 'Interpreter','latex');
grid on;

% 标签 A/B
text(1.1, 0, -0.1, 'A', 'FontSize', 20, 'Interpreter','latex', 'HorizontalAlignment','center');
text(-1.1, 0, 0, 'B', 'FontSize', 20, 'Interpreter','latex', 'HorizontalAlignment','center');

exportgraphics(gcf, 'PoincareSphere_3Trajectory_color_spline_2.8um_4.3um.png', 'Resolution',300)


%% 10 Test1 4.4um-5.9um L-16lines; Test2 6.0um-7.5um; Test3 13.5-15um

clc; clear;

% ---------------- 基本参数 ----------------
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.51078;

dn_dlambda = -0.0198418e6;        % /m
d2n_dlambda2 = 0.0492482e12;      % /m^2
v_phi = c / n;
k = omega0 / v_phi;  
k1 = (1/c) * (n - lambda0 * dn_dlambda); 
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2; 

t0 = 2e-14; 
Omega = 2 / t0; 
A0 = 1/2;
E0tilde = A0 / (Omega * sqrt(pi));
t3_values = linspace(-3*t0, 7*t0, 201);

% ---------------- Poincaré 球面设置 ----------------
figure;
[xs, ys, zs] = sphere(100);
colors4 = [1 0.6 0.6; 0.6 0.6 1; 0.6 1 0.6; 1 1 0.6];
C = zeros([size(ys), 3]);
for i = 1:size(ys,1)
    for j = 1:size(ys,2)
         y = ys(i,j); z = zs(i,j);
        if y >= 0 && z >= 0
            C(i,j,:) = colors4(1,:);
        elseif y < 0 && z >= 0
            C(i,j,:) = colors4(2,:);
        elseif y < 0 && z < 0
            C(i,j,:) = colors4(3,:);
        else
            C(i,j,:) = colors4(4,:);
        end
    end
end
ax = axes; hold(ax,'on');
surf(xs, ys, zs, C, 'EdgeColor','none', 'FaceAlpha',0.15, ...
     'HandleVisibility','off');

% ---- 轴线 ----
plot3([-1 1],[0 0],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, 'HandleVisibility','off');
plot3([0 0],[-1 1],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, 'HandleVisibility','off');
plot3([0 0],[0 0],[-1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, 'HandleVisibility','off');

% ---------------- 主循环：11 条轨迹 ----------------
L_list = linspace(13.5e-6, 15.0e-6, 16);
colors = jet(numel(L_list));
hLines = gobjects(numel(L_list),1);
legend_entries = cell(1,numel(L_list));


for li = 1:numel(L_list)
    L = L_list(li);
    Ex = A0 * exp(-t3_values.^2/t0^2) .* exp(-1i*omega0*t3_values) .* exp(1i*k0*L);
    Stokes = zeros(numel(t3_values),4);

    for idx = 1:numel(t3_values)
        t3 = t3_values(idx);
        A = 1/Omega^2 - 1i * k2 * L / 2;
        B = k1 * L - t3;
        Ey = E0tilde .* exp(1i*k*L) .* exp(-1i*omega0*t3) ...
            .* sqrt(pi./A) .* exp(-B.^2./(4*A));
        S0 = abs(Ex(idx)).^2 + abs(Ey).^2;
        S1 = abs(Ex(idx)).^2 - abs(Ey).^2;
        S2 = 2*real(Ex(idx).*conj(Ey));
        S3 = - 2*imag(Ex(idx).*conj(Ey));
        Stokes(idx,:) = [S1/S0, S2/S0, S3/S0, (S1^2+S2^2+S3^2)/S0^2];
    end

    % 插值平滑轨迹
    t_idx = 1:numel(t3_values);
    tt = linspace(1, numel(t3_values), 20*numel(t3_values));
    s1_spline = spline(t_idx, Stokes(:,1), tt);
    s2_spline = spline(t_idx, Stokes(:,2), tt);
    s3_spline = spline(t_idx, Stokes(:,3), tt);

    % 平滑轨迹线
    hLines(li) = plot3(s1_spline, s2_spline, s3_spline, ...
                       'Color', colors(li,:), 'LineWidth', 1.5);

    % 散点轨迹（不进 legend）
    scatter3(Stokes(:,1),Stokes(:,2),Stokes(:,3),20, ...
             linspace(0,1,numel(t3_values)), 'filled', ...
             'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.6, ...
             'HandleVisibility','off');

    legend_entries{li} = sprintf('$L = %.1f\\,\\mu\\mathrm{m}$', L*1e6);
    % --- 添加每条轨迹的中心点标注 ---
mid_idx = round(1.4*numel(t3_values)/3);
s1c = Stokes(mid_idx,1);
s2c = Stokes(mid_idx,2);
s3c = Stokes(mid_idx,3);
text(s1c, s2c, s3c, ...
     sprintf('%.1f\\,$\\mu$m', L*1e6), ...
     'FontSize', 10, 'Interpreter','latex', ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','middle', ...
     'Color', colors(li,:));

end

% ---------------- 图像设置 ----------------
view(120, 30);
%legend(hLines, legend_entries, 'Location','northeast', 'Box','on', 'Interpreter','latex');

axis equal;
set(gca, 'FontSize', 14);
xticks([-1, 0, 1]); yticks([-1, 0, 1]); zticks([-1, 0, 1]);
xlabel('$s_1$', 'FontSize', 24, 'Interpreter','latex');
ylabel('$s_2$', 'FontSize', 24, 'Interpreter','latex');
zlabel('$s_3$', 'FontSize', 24, 'Interpreter','latex');
grid on;

% 标签 A/B
text(1.1, 0, -0.1, 'A', 'FontSize', 20, 'Interpreter','latex', 'HorizontalAlignment','center');
text(-1.1, 0, 0, 'B', 'FontSize', 20, 'Interpreter','latex', 'HorizontalAlignment','center');

exportgraphics(gcf, 'PoincareSphere_3Trajectory_color_spline_13.5um_15um.png', 'Resolution',300)



%% 10Video L: 2.5um-20um
clc; clear all;

% ---------------- 基本参数 ----------------
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.51078;

dn_dlambda = -0.0198418e6;
d2n_dlambda2 = 0.0492482e12;
v_phi = c / n;
k = omega0 / v_phi;  
k1 = (1/c) * (n - lambda0 * dn_dlambda); 
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2; 

t0 = 2e-14; 
Omega = 2 / t0; 
A0 = 1/2;
E0tilde = A0 / (Omega * sqrt(pi));
t3_values = linspace(-1*t0, 3*t0, 41);

% ---------------- 创建视频 ----------------
videoFileName = 'poincare_trajectory_singlecolor.mp4';
v = VideoWriter(videoFileName, 'MPEG-4');
v.FrameRate = 5;
open(v);

% ---------------- Poincaré 球面设置 ----------------
figure('Color','w');  
%fig = figure('Color','w', 'Position', [100, 100, 1920/2, 1080/2]);  % Full HD 分辨率

[xs, ys, zs] = sphere(100);
colors4 = [1 0.6 0.6; 0.6 0.6 1; 0.6 1 0.6; 1 1 0.6];
C = zeros([size(ys), 3]);
for i = 1:size(ys,1)
    for j = 1:size(ys,2)
        y = ys(i,j); z = zs(i,j);
        if y >= 0 && z >= 0
            C(i,j,:) = colors4(1,:);
        elseif y < 0 && z >= 0
            C(i,j,:) = colors4(2,:);
        elseif y < 0 && z < 0
            C(i,j,:) = colors4(3,:);
        else
            C(i,j,:) = colors4(4,:);
        end
    end
end
ax = axes; hold(ax,'on');
surf(xs, ys, zs, C, 'EdgeColor','none', 'FaceAlpha',0.15);

plot3([-1 1],[0 0],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5);
plot3([0 0],[-1 1],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5);
plot3([0 0],[0 0],[-1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5);

view(120, 30);
axis equal;
set(gca, 'FontSize', 14);
xticks([-1, 0, 1]); yticks([-1, 0, 1]); zticks([-1, 0, 1]);
xlabel('$s_1$', 'FontSize', 24, 'Interpreter','latex');
ylabel('$s_2$', 'FontSize', 24, 'Interpreter','latex');
zlabel('$s_3$', 'FontSize', 24, 'Interpreter','latex');
grid on;
text(1.1, 0, -0.1, 'I', 'FontSize', 20, 'Interpreter','latex', 'HorizontalAlignment','center');
text(-1.1, 0.1, 0, 'F', 'FontSize', 20, 'Interpreter','latex', 'HorizontalAlignment','center');

% ---------------- 动画主循环 ----------------
L_list = linspace(2.5e-6, 20e-6, 200-25+1);
line_color = [0 0 1];      % 
scatter_color = [1 0 0];  % 灰色

for li = 1:numel(L_list)
    L = L_list(li);
    Ex = A0 * exp(-t3_values.^2/t0^2) .* exp(-1i*omega0*t3_values) .* exp(1i*k0*L);
    Stokes = zeros(numel(t3_values),4);

    for idx = 1:numel(t3_values)
        t3 = t3_values(idx);
        A = 1/Omega^2 - 1i * k2 * L / 2;
        B = k1 * L - t3;
        Ey = E0tilde .* exp(1i*k*L) .* exp(-1i*omega0*t3) ...  %k0？？？？？？？？？？
            .* sqrt(pi./A) .* exp(-B.^2./(4*A));
        S0 = abs(Ex(idx)).^2 + abs(Ey).^2;
        S1 = abs(Ex(idx)).^2 - abs(Ey).^2;
        S2 = 2*real(Ex(idx).*conj(Ey));
        S3 = - 2*imag(Ex(idx).*conj(Ey));
        Stokes(idx,:) = [S1/S0, S2/S0, S3/S0, (S1^2+S2^2+S3^2)/S0^2];
    end

    t_idx = 1:numel(t3_values);
    tt = linspace(1, numel(t3_values), 20*numel(t3_values));
    s1_spline = spline(t_idx, Stokes(:,1), tt);
    s2_spline = spline(t_idx, Stokes(:,2), tt);
    s3_spline = spline(t_idx, Stokes(:,3), tt);

    % 清除上一帧轨迹
    if exist('hPrev', 'var') && isvalid(hPrev)
        delete(hPrev);
    end
    if exist('hScatter', 'var')
        delete(hScatter);
    end

    % 绘制平滑轨迹线
    hPrev = plot3(s1_spline, s2_spline, s3_spline, ...
                  'Color', line_color, 'LineWidth', 2);

    % 绘制散点轨迹
    % hScatter = scatter3(Stokes(:,1), Stokes(:,2), Stokes(:,3), ...
    %                     20, scatter_color, 'filled', ...
    %                     'MarkerFaceAlpha', 0.6, ...
    %                     'MarkerEdgeColor', 'none');
   hScatter = scatter3(Stokes(:,1), Stokes(:,2), Stokes(:,3), ...
    30, linspace(0,1,numel(t3_values)), 'filled', ...
    'MarkerFaceAlpha', 0.9, ...
    'MarkerEdgeColor', 'none');
    colormap(jet); % 可选：设置颜色映射方案，比如 jet、parula 等
    


    % 更新标题
    title(sprintf('$L$ = %.2f $\\mu$m', L*1e6), ...
          'Interpreter', 'latex', 'FontSize', 20);
    
    %img = print(fig,'-RGBImage','-r150');  % 300 DPI 导出当前 figure 为图像
   %writeVideo(v, img);
  
    
    % 保存视频帧
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
disp(['✅ 视频已保存为: ', videoFileName]);
















%% 11 2D
% --------------------- 常量设置 ---------------------
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.5108;

% -------- 色散参数 --------
dn_dlambda = -0.0198e6;        % /m
d2n_dlambda2 = 0.0491e12;      % /m^2
v_phi = c / n;
k = omega0 / v_phi;
k1 = (1/c) * (n - lambda0 * dn_dlambda);
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

% -------- 脉冲参数 --------
t0 = 2e-14;
Omega = 2 / t0;
E0tilde = 1;
A0 = E0tilde * Omega * sqrt(pi);
t3_values = linspace(-1*t0, 3*t0, 41);

% ------------------ 图形设置 ------------------
figure; clf; hold on;


% ---- 绘制四色透明背景 ----
theta = linspace(0, pi/2, 100);
fill([0 cos(theta) 0], [0 sin(theta) 0], [1 0.6 0.6], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.15, 'HandleVisibility', 'off');  % 红右上

fill([0 -cos(theta) 0], [0 sin(theta) 0], [0.6 0.6 1], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.15, 'HandleVisibility', 'off');  % 蓝左上

fill([0 -cos(theta) 0], [0 -sin(theta) 0], [0.6 1 0.6], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.15, 'HandleVisibility', 'off');  % 绿左下

fill([0 cos(theta) 0], [0 -sin(theta) 0], [1 1 0.6], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.15, 'HandleVisibility', 'off');  % 黄右下

% ---- 单位圆轮廓 ----
%theta_full = linspace(0, 2*pi, 500);
%plot(cos(theta_full), sin(theta_full), 'k-', 'LineWidth', 1);

% ---- L 参数设置 ----
L_list = [5e-6, 10e-6, 15e-6];%5e-6, 10e-6, 15e-6
colors = {[0 1 1], [0 0 1], [1 0 0]};  % cyan, blue, red
legend_entries = cell(1, length(L_list));
plot_handles = gobjects(1, length(L_list));  % 保存句柄

for li = 1:length(L_list)
    L = L_list(li);
    Ex = A0 * exp(-t3_values.^2 / t0^2) .* exp(-1i * omega0 * t3_values) .* exp(1i * k0 * L);
    Stokes = zeros(length(t3_values), 5);
    
    for idx = 1:length(t3_values)
        t3 = t3_values(idx);
        A = 1/Omega^2 - 1i * k2 * L / 2;
        B = k1 * L - t3;
        Ey = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3) ...
            .* sqrt(pi ./ A) .* exp(-B.^2 ./ (4 * A));
        
        S0 = abs(Ex(idx)).^2 + abs(Ey).^2;
        S1 = abs(Ex(idx)).^2 - abs(Ey).^2;
        S2 = 2 * real(Ex(idx) .* conj(Ey));
        S3 =- 2 * imag(Ex(idx) .* conj(Ey));
        Stokes(idx, :) = [S1/S0, S2/S0, S3/S0, (S1^2+S2^2+S3^2)/S0^2, S2/S3];
    end

    % ---- 绘制轨迹线并保存句柄用于 legend ----
plot_handles(li) = plot(Stokes(:,2), Stokes(:,3), '-', ...
                        'Color', colors{li}, 'LineWidth', 2);

% ---- 添加透明散点但不加入 legend ----
scatter(Stokes(:,2), Stokes(:,3), ...
    20, linspace(0, 1, length(t3_values)), 'filled', ...
    'MarkerEdgeColor', 'none', ...
    'MarkerFaceAlpha', 0.6, ...
    'HandleVisibility', 'off');  % ✅ 不显示在 legend 中 'MarkerFaceColor', colors{li}, ...



    legend_entries{li} = sprintf('L = %.0f \\mum', L*1e6);
end

% ---------------- 图形美化 ----------------
axis equal;
xlim([-1.05, 1.05]);
ylim([-1.05, 1.05]);
set(gca, 'FontSize', 14);  %'FontName', 'Times New Roman',  
xlabel('$s_2$',  'FontSize', 24, 'Interpreter', 'latex');%
ylabel('$s_3$', 'FontSize', 24, 'Interpreter', 'latex');%
title('Polarization Trajectories on Poincaré Sphere (2D Projection)', 'FontWeight', 'bold');
legend(plot_handles, legend_entries, 'Location', 'northeast', 'Box', 'on');


% ---- 可选导出高分辨率图像 ----
% exportgraphics(gcf, 'Poincare2DProjection_3trace_color.png', 'Resolution', 300);




%% 20250901
%% 12 L=6400um; sphere - spline - values
% --------------------- Constants ---------------------
c = 3e8;
lambda0 = 800e-9;
omega0 = 2*pi*c / lambda0;
k0 = 2*pi / lambda0;
n = 1.51078;

% ---------------- Dispersion Parameters --------------
dn_dlambda = -0.0198418e6;        % /m
d2n_dlambda2 = 0.0492482e12;      % /m^2
v_phi = c / n;
k = omega0 / v_phi;  
k1 = (1/c) * (n - lambda0 * dn_dlambda); 
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2; 

% ------------------ Pulse Parameters ------------------
t0 = 2e-14; 
Omega = 2 / t0; 
%E0tilde = 1;
%A0 = E0tilde * Omega * sqrt(pi); 
A0 = 1/2;
E0tilde = A0/(Omega * sqrt(pi));
t3_values = linspace(3.250/2*10^3*t0, 3.270/2*10^3*t0, 101);

% ------------------ Plot Setup ------------------
figure;

% ---- Create colored unit sphere surface ----
[xs, ys, zs] = sphere(100);
colors4 = [1 0.6 0.6; 0.6 0.6 1; 0.6 1 0.6; 1 1 0.6];
C = zeros([size(ys), 3]);
for i = 1:size(ys,1)
    for j = 1:size(ys,2)
         y = ys(i,j); z = zs(i,j);
        if y >= 0 && z >= 0
            C(i,j,:) = colors4(1,:);
        elseif y < 0 && z >= 0
            C(i,j,:) = colors4(2,:);
        elseif y < 0 && z < 0
            C(i,j,:) = colors4(3,:);
        else
            C(i,j,:) = colors4(4,:);
        end
    end
end
ax = axes; hold(ax,'on');

% ---- Create colored unit sphere surface ----
surf(xs, ys, zs, C, 'EdgeColor','none', 'FaceAlpha',0.15, ...
     'HandleVisibility','off');  % 不进入 legend

% ---- Add axis lines ----
plot3([-1 1],[0 0],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, ...
      'HandleVisibility','off');
plot3([0 0],[-1 1],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, ...
      'HandleVisibility','off');
plot3([0 0],[0 0],[-1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5, ...
      'HandleVisibility','off');

% ---- Loop over different L values ----
L_list = 6.4e-3;

col = [0 1 1;   % cyan
       0 0 1;   % blue
       1 0 0;   % red
       0.6 0 0.6];% purple 

hLines = gobjects(numel(L_list),1);   % 保存三条线的句柄
legend_entries = cell(1,numel(L_list));

%t_delay_x = 32568*10^(-15); %t= 32568fs;
%t_delay_x = 32584*10^(-15);
%t_delay_x = 5088.85*6.4*10^(-15);
t_delay_x = 5088.85*6.4*10^(-15);

for li = 1:numel(L_list)
    L = L_list(li);
    Ex = A0 * exp(-(t3_values-t_delay_x).^2/t0^2) .* exp(-1i*omega0*t3_values) .* exp(1i*k0*L);
    
    Stokes = zeros(numel(t3_values),4);

    for idx = 1:numel(t3_values)
        t3 = t3_values(idx);
        A = 1/Omega^2 - 1i * k2 * L / 2;
        B = k1 * L - t3;
         Bt(idx) = k1 * L - t3;
        Ey = E0tilde .* exp(1i*k*L) .* exp(-1i*omega0*t3) ...
            .* sqrt(pi./A) .* exp(-B.^2./(4*A));
        S0 = abs(Ex(idx)).^2 + abs(Ey).^2;
        S1 = abs(Ex(idx)).^2 - abs(Ey).^2;
        S2 = 2*real(Ex(idx).*conj(Ey));
        S3 = - 2*imag(Ex(idx).*conj(Ey));
        Stokes(idx,:) = [S1/S0, S2/S0, S3/S0, (S1^2+S2^2+S3^2)/S0^2];
    end

    % 使用样条插值进行平滑
    t_idx = 1:numel(t3_values);                     % 原始点参数
    tt = linspace(1, numel(t3_values), 20*numel(t3_values));  % 插值参数

    % 对每个分量进行样条插值
    s1_spline = spline(t_idx, Stokes(:,1), tt);
    s2_spline = spline(t_idx, Stokes(:,2), tt);
    s3_spline = spline(t_idx, Stokes(:,3), tt);

    % 平滑轨迹实线（进入 legend）
    hLines(li) = plot3(s1_spline, s2_spline, s3_spline, ...
                   'Color', col(li,:), 'LineWidth', 2);


    % 彩色散点（不进入 legend）
    scatter3(Stokes(:,1),Stokes(:,2),Stokes(:,3), 30, ...
             linspace(0,1,numel(t3_values)), 'filled', ...
             'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.9, ...
             'HandleVisibility','off');

    % % 彩色散点（颜色与线条匹配，不进入 legend）
    %      scatter3(Stokes(:,1), Stokes(:,2), Stokes(:,3), 20, ...
    %      repmat(col(li,:), size(Stokes,1), 1), 'filled', ...
    %      'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.6, ...
    %      'HandleVisibility','off');


    legend_entries{li} = sprintf('$L =$ %.0f $\\mu\\mathrm{m}$', L*1e6);
end
%view(120, 45);
view(120, 30);
%view(0, 90);%from top view
%legend(hLines, legend_entries, 'Location','northeast', 'Box', 'on', 'Interpreter', 'latex');

% ---- Final Plot Settings ----
axis equal;
set(gca, 'FontSize', 14);  %'FontName', 'Times New Roman',  
xticks([-1, 0, 1]);
yticks([-1, 0, 1]);
zticks([-1, 0, 1]);
xlabel('$s_1$', 'FontSize', 24, 'Interpreter', 'latex');
ylabel('$s_2$', 'FontSize', 24, 'Interpreter', 'latex');
zlabel('$s_3$', 'FontSize', 24, 'Interpreter', 'latex');
%title('Polarization Trajectories on Poincaré Sphere for Different L');
grid on;

text(1.1, 0, -0.1, 'I', 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
text(-1.1, 0.1, 0, 'F', 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');


exportgraphics(gcf, 'PoincareSphere_return_spline.png', 'Resolution',300)


%% 13  L=6400um; sphere - spline - t_delay 20 values
% ------------------ Pulse Parameters ------------------
t0 = 2e-14; 
Omega = 2 / t0; 
A0 = 1/2;
E0tilde = A0/(Omega * sqrt(pi));
t3_values = linspace(3.250/2*10^3*t0, 3.270/2*10^3*t0, 101);

% ------------------ Plot Setup ------------------
figure;

% ---- Create colored unit sphere surface ----
[xs, ys, zs] = sphere(100);
colors4 = [1 0.6 0.6; 0.6 0.6 1; 0.6 1 0.6; 1 1 0.6];
C = zeros([size(ys), 3]);
for i = 1:size(ys,1)
    for j = 1:size(ys,2)
        y = ys(i,j); z = zs(i,j);
        if y >= 0 && z >= 0
            C(i,j,:) = colors4(1,:);
        elseif y < 0 && z >= 0
            C(i,j,:) = colors4(2,:);
        elseif y < 0 && z < 0
            C(i,j,:) = colors4(3,:);
        else
            C(i,j,:) = colors4(4,:);
        end
    end
end
ax = axes; hold(ax,'on');

surf(xs, ys, zs, C, 'EdgeColor','none', 'FaceAlpha',0.15, 'HandleVisibility','off');

plot3([-1 1],[0 0],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5,'HandleVisibility','off');
plot3([0 0],[-1 1],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5,'HandleVisibility','off');
plot3([0 0],[0 0],[-1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5,'HandleVisibility','off');

% ---- Loop over different L values ----
L_list = 6.4e-3;

% ===================== 修改处（开始）=====================
% 生成 32490 到 32650（fs）的 20 个等差点，并换算为秒
t_delay_list = linspace(32490, 32650, 20) * 1e-15;   % [s]

% 为不同 t_delay 配色（20 色）
cols = lines(numel(t_delay_list));

% 预分配图例（如果只想显示 L 的图例也可以简化）
hLines = gobjects(numel(L_list)*numel(t_delay_list),1);
legend_entries = cell(1, numel(L_list)*numel(t_delay_list));
% ===================== 修改处（结束）=====================

for li = 1:numel(L_list)
    L = L_list(li);

    % 对每个 t_delay_x 循环
    for mi = 1:numel(t_delay_list)
        % ===================== 修改处：把原来的单值改为列表取值 ======================
        t_delay_x = t_delay_list(mi);   % 当前延时（秒）
        % =========================================================================

        Ex = A0 * exp(-(t3_values - t_delay_x).^2/t0^2) .* exp(-1i*omega0*t3_values) .* exp(1i*k0*L);

        Stokes = zeros(numel(t3_values),4);

        for idx = 1:numel(t3_values)
            t3 = t3_values(idx);
            A = 1/Omega^2 - 1i * k2 * L / 2;
            B = k1 * L - t3;
            Ey = E0tilde .* exp(1i*k*L) .* exp(-1i*omega0*t3) ...
                .* sqrt(pi./A) .* exp(-B.^2./(4*A));

            S0 = abs(Ex(idx)).^2 + abs(Ey).^2;
            S1 = abs(Ex(idx)).^2 - abs(Ey).^2;
            S2 = 2*real(Ex(idx).*conj(Ey));
            S3 = - 2*imag(Ex(idx).*conj(Ey));
            Stokes(idx,:) = [S1/S0, S2/S0, S3/S0, (S1^2+S2^2+S3^2)/S0^2];
        end

        % 样条平滑
        t_idx = 1:numel(t3_values);
        tt = linspace(1, numel(t3_values), 20*numel(t3_values));
        s1_spline = spline(t_idx, Stokes(:,1), tt);
        s2_spline = spline(t_idx, Stokes(:,2), tt);
        s3_spline = spline(t_idx, Stokes(:,3), tt);

        % 绘制
        lineIdx = (li-1)*numel(t_delay_list) + mi;
        hLines(lineIdx) = plot3(s1_spline, s2_spline, s3_spline, ...
                                'Color', cols(mi,:), 'LineWidth', 2);

        scatter3(Stokes(:,1),Stokes(:,2),Stokes(:,3), 30, ...
                 linspace(0,1,numel(t3_values)), 'filled', ...
                 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.9, ...
                 'HandleVisibility','off');

        % 图例（包含 t_delay，单位 fs）
        legend_entries{lineIdx} = sprintf('$L=%.0f\\,\\mu\\mathrm{m}$, $t_d=%.2f\\,\\mathrm{fs}$', ...
                                          L*1e6, t_delay_x*1e15);
    end
end

view(120, 30);
legend(hLines, legend_entries, 'Location','northeast', 'Box','on', 'Interpreter','latex');

axis equal;
set(gca, 'FontSize', 14);
xticks([-1, 0, 1]); yticks([-1, 0, 1]); zticks([-1, 0, 1]);
xlabel('$s_1$', 'FontSize', 24, 'Interpreter', 'latex');
ylabel('$s_2$', 'FontSize', 24, 'Interpreter', 'latex');
zlabel('$s_3$', 'FontSize', 24, 'Interpreter', 'latex');
grid on;

text(1.1, 0, -0.1, 'A', 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
text(-1.1, 0.1, 0, 'B', 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');


exportgraphics(gcf, 'PoincareSphere_return_spline.png', 'Resolution',300)


%% 14 Simple Video L=6400um; sphere - spline - t_delay 161 values
% ------------------ 视频设置 ------------------
v = VideoWriter('poincare_trajectory.mp4','MPEG-4');  % 输出 MP4 文件
v.FrameRate = 5;   % 每秒 5 帧，可根据需要调整
open(v);

% ------------------ Pulse Parameters ------------------
t0 = 2e-14; 
Omega = 2 / t0; 
A0 = 1/2;
E0tilde = A0/(Omega * sqrt(pi));
t3_values = linspace(3.250/2*10^3*t0, 3.270/2*10^3*t0, 41);

% ---- 固定参数 ----
L = 6.4e-3;

% 延时范围（20 个点，fs -> s）
t_delay_list = linspace(32490, 32650, 161) * 1e-15;   % [s]

% ------------------ 循环绘制并保存 ------------------
for mi = 1:numel(t_delay_list)
    t_delay_x = t_delay_list(mi);

    Ex = A0 * exp(-(t3_values - t_delay_x).^2/t0^2) .* ...
         exp(-1i*omega0*t3_values) .* exp(1i*k0*L);

    Stokes = zeros(numel(t3_values),4);
    for idx = 1:numel(t3_values)
        t3 = t3_values(idx);
        A = 1/Omega^2 - 1i * k2 * L / 2;
        B = k1 * L - t3;
        Ey = E0tilde .* exp(1i*k*L) .* exp(-1i*omega0*t3) ...
            .* sqrt(pi./A) .* exp(-B.^2./(4*A));

        S0 = abs(Ex(idx)).^2 + abs(Ey).^2;
        S1 = abs(Ex(idx)).^2 - abs(Ey).^2;
        S2 = 2*real(Ex(idx).*conj(Ey));
        S3 = - 2*imag(Ex(idx).*conj(Ey));
        Stokes(idx,:) = [S1/S0, S2/S0, S3/S0, (S1^2+S2^2+S3^2)/S0^2];
    end

    % 样条平滑
    t_idx = 1:numel(t3_values);
    tt = linspace(1, numel(t3_values), 20*numel(t3_values));
    s1_spline = spline(t_idx, Stokes(:,1), tt);
    s2_spline = spline(t_idx, Stokes(:,2), tt);
    s3_spline = spline(t_idx, Stokes(:,3), tt);

    % ------------------ 绘制 ------------------
    clf;   % 清空图
    hold on;
    axis equal;
    grid on;

    % Poincaré sphere
    [xs, ys, zs] = sphere(100);
    surf(xs, ys, zs, 'FaceAlpha',0.05, 'EdgeColor','none');

    % 坐标轴
    plot3([-1 1],[0 0],[0 0],'k--');
    plot3([0 0],[-1 1],[0 0],'k--');
    plot3([0 0],[0 0],[-1 1],'k--');

    % 曲线
    plot3(s1_spline, s2_spline, s3_spline, 'r', 'LineWidth', 2);
    scatter3(Stokes(:,1),Stokes(:,2),Stokes(:,3), 30, ...
             linspace(0,1,numel(t3_values)), 'filled');

    view(120,30);
    xlabel('$s_1$', 'Interpreter','latex');
    ylabel('$s_2$', 'Interpreter','latex');
    zlabel('$s_3$', 'Interpreter','latex');
    title(sprintf('t_{delay} = %.1f fs', t_delay_x*1e15));

    % ------------------ 保存帧 ------------------
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% ------------------ 关闭视频 ------------------
close(v);




%% 14 GGood video
% ===================== 常量与色散参数 =====================
c = 3e8;
lambda0 = 800e-9;
omega0  = 2*pi*c / lambda0;
k0      = 2*pi / lambda0;
n       = 1.51078;

% ---------------- Dispersion Parameters --------------
dn_dlambda   = -0.0198418e6;       % /m
d2n_dlambda2 = 0.0492482e12;       % /m^2
v_phi = c / n;
k     = omega0 / v_phi;
k1    = (1/c) * (n - lambda0 * dn_dlambda);
k2    = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

% ===================== 脉冲与扫描参数 =====================
t0        = 2e-14;
Omega     = 2 / t0;
A0        = 1/2;
E0tilde   = A0/(Omega * sqrt(pi));
t3_values = linspace(3.250/2*10^3*t0, 3.270/2*10^3*t0, 41);
N         = numel(t3_values);

% 传播长度
L = 6.4e-3;

% 扫描延时（fs -> s）
t_delay_list = linspace(32490, 32650, 161) * 1e-15;

% ===================== 预生成四象限配色庞加莱球 =====================
[xs, ys, zs] = sphere(100);
colors4 = [1   0.6 0.6;   % Q1 (y>=0, z>=0)
           0.6 0.6 1;     % Q2 (y<0,  z>=0)
           0.6 1   0.6;   % Q3 (y<0,  z<0)
           1   1   0.6];  % Q4 (y>=0, z<0)
C = zeros([size(ys), 3]);
for i = 1:size(ys,1)
    for j = 1:size(ys,2)
        y = ys(i,j); z = zs(i,j);
        if y >= 0 && z >= 0
            C(i,j,:) = colors4(1,:);
        elseif y < 0 && z >= 0
            C(i,j,:) = colors4(2,:);
        elseif y < 0 && z < 0
            C(i,j,:) = colors4(3,:);
        else
            C(i,j,:) = colors4(4,:);
        end
    end
end

% ===================== 视频设置 =====================
v = VideoWriter('poincare_trajectory_2.mp4','MPEG-4');
v.FrameRate = 5;
v.Quality   = 100;
open(v);

line_color = [0 0 1];      % 轨迹主色

% ===================== 图窗与静态元素（只建一次） =====================
fig = figure('Color','w');%,'Renderer','opengl'
ax  = axes('Parent',fig); hold(ax,'on');
view(ax,120,30);
axis(ax,'equal'); axis(ax,[-1 1 -1 1 -1 1.08]); grid(ax,'on');
set(ax, 'FontSize', 14);
xticks(ax,[-1 0 1]); yticks(ax,[-1 0 1]); zticks(ax,[-1 0 1]);

% 彩色四象限庞加莱球
surf(ax, xs, ys, zs, C, 'EdgeColor','none', 'FaceAlpha',0.15, 'HandleVisibility','off');

% 灰色虚线坐标轴
plot3(ax,[-1 1],[0 0],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5,'HandleVisibility','off');
plot3(ax,[0 0],[-1 1],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5,'HandleVisibility','off');
plot3(ax,[0 0],[0 0],[-1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5,'HandleVisibility','off');

xlabel(ax,'$s_1$','FontSize',24,'Interpreter','latex');
ylabel(ax,'$s_2$','FontSize',24,'Interpreter','latex');
zlabel(ax,'$s_3$','FontSize',24,'Interpreter','latex');

% A / B 标注
text(ax,  1.1, 0,  -0.1, 'A', 'FontSize',20,'Interpreter','latex','HorizontalAlignment','center');
text(ax, -1.1, 0.1, 0,    'B', 'FontSize',20,'Interpreter','latex','HorizontalAlignment','center');

% ===== 散点与轨迹图元：只创建一次 =====
hLine = plot3(ax, nan, nan, nan, 'Color', line_color, 'LineWidth', 2);

hScatter = scatter3(ax, nan, nan, nan, 30, nan, ...
    'filled', 'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor', 'none', ...
    'HandleVisibility','off');

% 颜色映射：使用 jet，并固定 [0,1]，加 colorbar（可选）
colormap(ax,'jet');
caxis(ax,[0 1]);
%colorbar('eastoutside');

% 标题句柄
ht = title(ax,'','Interpreter','latex','FontSize',24);

% ===================== 预计算（向量化，避免每帧重复） =====================
A        = 1/Omega^2 - 1i * k2 * L / 2;                     % 常量
constEy  = E0tilde .* exp(1i*k*L) .* sqrt(pi./A);           % 常量（含 sqrt）
phase_t3 = exp(-1i*omega0*t3_values);                       % e^{-i ω0 t3}
B_vec    = k1 * L - t3_values;                              % 向量
Ey_base  = constEy .* phase_t3 .* exp(-(B_vec.^2)./(4*A));  % 与 t_delay 无关

% 样条平滑索引
tt    = linspace(1, N, 100*N);
t_idx = 1:N;

% ===================== 帧循环（只更新数据） =====================
for mi = 1:numel(t_delay_list)
    t_delay_x = t_delay_list(mi);

    % ------ 场（向量化） ------
    Ex = A0 .* exp(-((t3_values - t_delay_x).^2)/t0^2) .* phase_t3 .* exp(1i*k0*L);
    Ey = Ey_base;  % 与 t_delay 无关（按你的解析式）

    % ------ 斯托克斯（向量化 + 安全归一） ------
    Ex2 = abs(Ex).^2; Ey2 = abs(Ey).^2;
    S0  = max(Ex2 + Ey2, eps);      % 防止除零
    S1  = Ex2 - Ey2;
    S2  = 2*real(Ex .* conj(Ey));
    S3  = -2*imag(Ex .* conj(Ey));

    % 归一化后的 s1,s2,s3
    s1n = S1 ./ S0; s2n = S2 ./ S0; s3n = S3 ./ S0;

    % 组装成 Stokes(.,1:3) 便于与你的写法一致
    Stokes = [s1n(:), s2n(:), s3n(:)];

    % ------ 样条平滑（用于轨迹线） ------
    s1_spline = spline(t_idx, s1n, tt);
    s2_spline = spline(t_idx, s2n, tt);
    s3_spline = spline(t_idx, s3n, tt);

    % ------ 更新轨迹线 ------
    set(hLine, 'XData', s1_spline, 'YData', s2_spline, 'ZData', s3_spline, ...
               'Color', line_color, 'LineWidth', 2);

    % ------ 更新散点（使用 jet 映射到 [0,1]） ------
    set(hScatter, ...
        'XData', Stokes(:,1), ...
        'YData', Stokes(:,2), ...
        'ZData', Stokes(:,3), ...
        'CData', linspace(0,1,N), ...  % 颜色随样本索引渐变
        'SizeData', 30, ...
        'MarkerFaceAlpha', 0.9, ...
        'MarkerEdgeColor', 'none');

    % 标题：当前延时
    set(ht, 'String', sprintf('$\\Delta t = %.1f\\,\\mathrm{fs}$', t_delay_x*1e15));

    drawnow limitrate;
  img = print(fig,'-RGBImage','-r100');  % 300 DPI 导出当前 figure 为图像
   writeVideo(v, img);
    % ------ 写入视频帧 ------
    %frame = getframe(fig);
    %writeVideo(v, frame);
end

% ===================== 收尾 =====================
close(v);
disp('✅ 已生成视频：poincare_trajectory_2.mp4');


%% 15 Combined: Ex/Ey (time trace) + Poincaré trajectory in one video
% ---------- Constants ----------
c = 3e8;
lambda0 = 800e-9;
omega0  = 2*pi*c / lambda0;
k0      = 2*pi / lambda0;
n       = 1.51078;

% ---------- Dispersion ----------
dn_dlambda   = -0.0198418e6;       % /m
d2n_dlambda2 = 0.0492482e12;       % /m^2
v_phi = c / n;
k     = omega0 / v_phi;
k1    = (1/c) * (n - lambda0 * dn_dlambda);
k2    = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

% ---------- Pulse / scan params ----------
t0      = 2e-14;         % 20 fs
Omega   = 2 / t0;
A0      = 1/2;
E0tilde = A0/(Omega * sqrt(pi));
L       = 6.4e-3;        % propagation distance

% Delay list (fs -> s)
t_delay_list = linspace(32490, 32650, 161) * 1e-15;

% ---------- Time arrays ----------
t3_time   = linspace(3.24/2*10^3*t0, 3.275/2*10^3*t0, 5000);   % 左图
t3_values = linspace(3.250/2*10^3*t0, 3.270/2*10^3*t0, 41);    % 右图
N         = numel(t3_values);

% ---------- Precompute common terms ----------
A_left   = 1/Omega^2 - 1i * k2 * L / 2;
B_left   = k1 * L - t3_time;
Ey_left0 = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3_time) ...
           .* sqrt(pi ./ A_left) .* exp(-B_left.^2 ./ (4*A_left));
Ey_left  = 2 * real(Ey_left0);

A_right   = 1/Omega^2 - 1i * k2 * L / 2;
constEy   = E0tilde .* exp(1i*k*L) .* sqrt(pi./A_right);
phase_t3  = exp(-1i*omega0*t3_values);
B_vec     = k1 * L - t3_values;
Ey_base   = constEy .* phase_t3 .* exp(-(B_vec.^2)./(4*A_right));

tt    = linspace(1, N, 100*N);
t_idx = 1:N;

% ---------- Colors ----------
custom_colors = [
    0.0, 0.9, 0.9;
    0,   0,   1;
    1,   0,   0
];
line_color = [0 0 1];

% Combined: Ex/Ey (time trace) + Poincaré trajectory in one video
% ...（前面常量和预计算部分保持不变）...

% ---------- Video ----------
v = VideoWriter('combined_pulse_poincare2.mp4','MPEG-4');
v.FrameRate = 20;
v.Quality   = 100;
open(v);

% ---------- Figure & layout ----------
fig = figure('Color','w','Position',[50 50 1300 600]);

% 用 1x6 网格：左边 2 格，右边 4 格
tiled = tiledlayout(fig,1,9,'Padding','compact','TileSpacing','compact');

% ===== Left panel: Ex/Ey (小) =====
ax1 = nexttile(tiled,[1 4]); hold(ax1,'on'); grid(ax1,'on');
set(ax1,'FontSize',14);
xlabel(ax1,'$t$ [fs]','Interpreter','latex','FontSize',16);
ylabel(ax1,'$E(t)$','Interpreter','latex','FontSize',16);

hEx = plot(ax1, t3_time*1e15, nan(size(t3_time)), 'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
hEy = plot(ax1, t3_time*1e15, Ey_left,            'Color',custom_colors(2,:), 'LineWidth',1.2);
legend(ax1, {'$E_x$','$E_y$'}, 'Interpreter','latex','FontSize',14,'Location','best');

% ===== Right panel: Poincaré sphere (大) =====
ax2 = nexttile(tiled,[1 5]); hold(ax2,'on');
view(ax2,120,30); axis(ax2,'equal'); axis(ax2,[-1 1 -1 1 -1 1.08]); grid(ax2,'on');
set(ax2,'FontSize',14);
xticks(ax2,[-1 0 1]); yticks(ax2,[-1 0 1]); zticks(ax2,[-1 0 1]);
xlabel(ax2,'$s_1$','FontSize',24,'Interpreter','latex');
ylabel(ax2,'$s_2$','FontSize',24,'Interpreter','latex');
zlabel(ax2,'$s_3$','FontSize',24,'Interpreter','latex');

% 彩色四象限球
[xs, ys, zs] = sphere(100);
C = zeros([size(ys),3]);
for i = 1:size(ys,1)
    for j = 1:size(ys,2)
        y = ys(i,j); z = zs(i,j);
        if y >= 0 && z >= 0
            C(i,j,:) = [1 0.6 0.6];
        elseif y < 0 && z >= 0
            C(i,j,:) = [0.6 0.6 1];
        elseif y < 0 && z < 0
            C(i,j,:) = [0.6 1 0.6];
        else
            C(i,j,:) = [1 1 0.6];
        end
    end
end
surf(ax2, xs, ys, zs, C, 'EdgeColor','none', 'FaceAlpha',0.15);

% 灰虚线坐标轴
plot3(ax2,[-1 1],[0 0],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5);
plot3(ax2,[0 0],[-1 1],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5);
plot3(ax2,[0 0],[0 0],[-1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5);

% 标记 A/B
text(ax2,  1.1, 0,  -0.1, 'I', 'FontSize',20,'Interpreter','latex','HorizontalAlignment','center');
text(ax2, -1.1, 0.1, 0,    'F', 'FontSize',20,'Interpreter','latex','HorizontalAlignment','center');

hLine    = plot3(ax2, nan, nan, nan, 'Color', line_color, 'LineWidth', 2);
hScatter = scatter3(ax2, nan, nan, nan, 30, nan, 'filled', ...
    'MarkerFaceAlpha',0.9, 'MarkerEdgeColor','none');

colormap(ax2,'jet'); caxis(ax2,[0 1]);

% ======= 全局标题（一个就够） =======
ht = title(tiled,'','Interpreter','latex','FontSize',20);

% ---------- Frame loop ----------
for mi = 1:numel(t_delay_list)
    t_delay_x = t_delay_list(mi);

    % 左边
    Ex_left_c = A0 * exp(-(t3_time - t_delay_x).^2 / t0^2) .* exp(-1i*omega0*t3_time);
    Ex_left   = 2 * real(Ex_left_c);
    set(hEx,'YData',Ex_left);
    set(hEy,'YData',Ey_left);

    % 右边
    Ex_r = A0 .* exp(-((t3_values - t_delay_x).^2)/t0^2) .* exp(-1i*omega0*t3_values) .* exp(1i*k0*L);
    Ey_r = Ey_base;
    Ex2 = abs(Ex_r).^2; Ey2 = abs(Ey_r).^2;
    S0  = max(Ex2 + Ey2, eps);
    S1  = Ex2 - Ey2;
    S2  = 2*real(Ex_r .* conj(Ey_r));
    S3  = -2*imag(Ex_r .* conj(Ey_r));
    s1n = S1 ./ S0; s2n = S2 ./ S0; s3n = S3 ./ S0;
    s1_spline = spline(t_idx, s1n, tt);
    s2_spline = spline(t_idx, s2n, tt);
    s3_spline = spline(t_idx, s3n, tt);
    set(hLine,'XData',s1_spline,'YData',s2_spline,'ZData',s3_spline);
    set(hScatter,'XData',s1n,'YData',s2n,'ZData',s3n,'CData',linspace(0,1,N));

    % 标题
    set(ht,'String',sprintf('$\\Delta t = %.1f\\,\\mathrm{fs}$', t_delay_x*1e15));

    drawnow limitrate;
    img = print(fig,'-RGBImage','-r100');  % 高分辨率帧
    writeVideo(v, img);
end

close(v);
disp('Saved video: combined_pulse_poincare.mp4');

%% 15 Good Combined: Ex/Ey (time trace) + Poincaré trajectory in one video (4:5 layout, left 2:3 aspect)
% ---------- Constants ----------
c = 3e8;
lambda0 = 800e-9;
omega0  = 2*pi*c / lambda0;
k0      = 2*pi / lambda0;
n       = 1.51078;

% ---------- Dispersion ----------
dn_dlambda   = -0.0198418e6;       % /m
d2n_dlambda2 = 0.0492482e12;       % /m^2
v_phi = c / n;
k     = omega0 / v_phi;
k1    = (1/c) * (n - lambda0 * dn_dlambda);
k2    = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

% ---------- Pulse / scan params ----------
t0      = 2e-14;         % 20 fs
Omega   = 2 / t0;
A0      = 1/2;
E0tilde = A0/(Omega * sqrt(pi));
L       = 6.4e-3;        % propagation distance

% Delay list (fs -> s)
t_delay_list = linspace(32490, 32650, 161) * 1e-15;

% ---------- Time arrays ----------
% Left panel waveform (dense)
t3_time   = linspace(3.24/2*10^3*t0, 3.275/2*10^3*t0, 5000);
% Right panel Stokes sampling (sparser)
t3_values = linspace(3.250/2*10^3*t0, 3.270/2*10^3*t0, 101);
N         = numel(t3_values);

% ---------- Precompute common terms ----------
% Left Ey (independent of delay)
A_left   = 1/Omega^2 - 1i * k2 * L / 2;
B_left   = k1 * L - t3_time;
Ey_left0 = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3_time) ...
           .* sqrt(pi ./ A_left) .* exp(-B_left.^2 ./ (4*A_left));
Ey_left  = 2 * real(Ey_left0);

% Right Ey base (independent of delay)
A_right   = 1/Omega^2 - 1i * k2 * L / 2;
constEy   = E0tilde .* exp(1i*k*L) .* sqrt(pi./A_right);
phase_t3  = exp(-1i*omega0*t3_values);
B_vec     = k1 * L - t3_values;
Ey_base   = constEy .* phase_t3 .* exp(-(B_vec.^2)./(4*A_right));

% Spline param for smooth trajectory
tt    = linspace(1, N, 20*N);
t_idx = 1:N;

% ---------- Colors ----------
custom_colors = [
    0.0, 0.9, 0.9;  % cyan
    0,   0,   1;    % blue
    1,   0,   0     % red
];
line_color = [0 0 1];

% ---------- Video ----------
v = VideoWriter('combined_pulse_poincare5.mp4','MPEG-4');
v.FrameRate = 5;
v.Quality   = 100;
open(v);

% ---------- Figure & layout ----------
fig = figure('Color','w','Position',[50 50 1350 600]);

% 1x9 网格：左 4，右 5 → 总体 4:5 比例
tiled = tiledlayout(fig,1,9,'Padding','compact','TileSpacing','compact');

% ===== Left panel: Ex/Ey (4/9 宽，纵横比 2:3) =====
ax1 = nexttile(tiled,[1 4]); hold(ax1,'on'); grid(ax1,'on');
set(ax1,'FontSize',14);
xlabel(ax1,'Time [fs]','Interpreter','latex','FontSize',16);
ylabel(ax1,'$E(t)$','Interpreter','latex','FontSize',16);
% 强制 2:3 的纵横比（宽:高 = 3:2）
pbaspect(ax1,[3 2 1]);

hEx = plot(ax1, t3_time*1e15, nan(size(t3_time)), 'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
hEy = plot(ax1, t3_time*1e15, Ey_left,            'Color',custom_colors(2,:), 'LineWidth',1.2);
legend(ax1, {'$E_x$','$E_y$'}, 'Interpreter','latex','FontSize',14,'Location','best');

% ===== Right panel: Poincaré sphere (5/9 宽) =====
ax2 = nexttile(tiled,[1 5]); hold(ax2,'on');
view(ax2,120,30);
axis(ax2,'equal'); axis(ax2,[-1 1 -1 1 -1 1]); grid(ax2,'on');
set(ax2,'FontSize',14);
xticks(ax2,[-1 0 1]); yticks(ax2,[-1 0 1]); zticks(ax2,[-1 0 1]);
xlabel(ax2,'$s_1$','FontSize',24,'Interpreter','latex');
ylabel(ax2,'$s_2$','FontSize',24,'Interpreter','latex');
zlabel(ax2,'$s_3$','FontSize',24,'Interpreter','latex');

% 彩色四象限球
[xs, ys, zs] = sphere(100);
colors4 = [1   0.6 0.6;   % Q1 (y>=0, z>=0)
           0.6 0.6 1;     % Q2 (y<0,  z>=0)
           0.6 1   0.6;   % Q3 (y<0,  z<0)
           1   1   0.6];  % Q4 (y>=0, z<0)
C = zeros([size(ys),3]);
for i = 1:size(ys,1)
    for j = 1:size(ys,2)
        y = ys(i,j); z = zs(i,j);
        if y >= 0 && z >= 0
            C(i,j,:) = colors4(1,:);
        elseif y < 0 && z >= 0
            C(i,j,:) = colors4(2,:);
        elseif y < 0 && z < 0
            C(i,j,:) = colors4(3,:);
        else
            C(i,j,:) = colors4(4,:);
        end
    end
end
surf(ax2, xs, ys, zs, C, 'EdgeColor','none', 'FaceAlpha',0.15);

% 灰虚线坐标轴与标记
plot3(ax2,[-1 1],[0 0],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5);
plot3(ax2,[0 0],[-1 1],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5);
plot3(ax2,[0 0],[0 0],[-1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',0.5);
text(ax2,  1.1, 0,  -0.1, 'I', 'FontSize',20,'Interpreter','latex','HorizontalAlignment','center');
text(ax2, -1.1, 0.1, 0,   'F', 'FontSize',20,'Interpreter','latex','HorizontalAlignment','center');

% 轨迹与散点
hLine    = plot3(ax2, nan, nan, nan, 'Color', line_color, 'LineWidth', 2);
hScatter = scatter3(ax2, nan, nan, nan, 35, nan, 'filled', ...
    'MarkerFaceAlpha',0.8, 'MarkerEdgeColor','none');
colormap(ax2,'jet'); clim(ax2,[0 1]); %caxis

% ========= 全局标题（一个） =========
ht = title(tiled,'','Interpreter','latex','FontSize',20);

% ---------- Frame loop ----------
for mi = 1:numel(t_delay_list)
    t_delay_x = t_delay_list(mi);

    % ---- Left panel fields ----
    Ex_left_c = A0 * exp(-(t3_time - t_delay_x).^2 / t0^2) .* exp(-1i*omega0*t3_time);
    Ex_left   = 2 * real(Ex_left_c);
    set(hEx,'YData',Ex_left);
    set(hEy,'YData',Ey_left);

    % ---- Right panel Stokes ----
    Ex_r = A0 .* exp(-((t3_values - t_delay_x).^2)/t0^2) .* exp(-1i*omega0*t3_values) .* exp(1i*k0*L);
    Ey_r = Ey_base;

    Ex2 = abs(Ex_r).^2; 
    Ey2 = abs(Ey_r).^2;
    S0  = max(Ex2 + Ey2, eps);
    S1  = Ex2 - Ey2;
    S2  = 2*real(Ex_r .* conj(Ey_r));
    S3  = -2*imag(Ex_r .* conj(Ey_r));
    s1n = S1 ./ S0; s2n = S2 ./ S0; s3n = S3 ./ S0;

    % Smooth trajectory
    s1_spline = spline(t_idx, s1n, tt);
    s2_spline = spline(t_idx, s2n, tt);
    s3_spline = spline(t_idx, s3n, tt);

    set(hLine,'XData',s1_spline,'YData',s2_spline,'ZData',s3_spline);
    set(hScatter,'XData',s1n,'YData',s2n,'ZData',s3n,'CData',linspace(0,1,N));

    % ---- Update global title ----
    set(ht,'String',sprintf('$\\Delta t = %.1f\\,\\mathrm{fs}$', t_delay_x*1e15));

    drawnow limitrate;
    img = print(fig,'-RGBImage','-r100');   % 高分辨率导出当前帧
    writeVideo(v, img);
end

% ---------- Finish ----------
close(v);
disp('Saved video: combined_pulse_poincare.mp4');

%% 16  Good!! Polarization: pulse waveform (left) + Poincaré trajectory (right)
% ---------- Constants ----------
c = 3e8;
lambda0 = 800e-9;
omega0  = 2*pi*c / lambda0;
k0      = 2*pi / lambda0;
n       = 1.51078;

% ---------- Dispersion ----------
dn_dlambda   = -0.0198418e6;       % /m
d2n_dlambda2 = 0.0492482e12;       % /m^2
v_phi = c / n;
k     = omega0 / v_phi;
k1    = (1/c) * (n - lambda0 * dn_dlambda);
k2    = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

% ---------- Pulse / scan params ----------
t0      = 2e-14;         % 20 fs
Omega   = 2 / t0;
A0      = 1/2;
E0tilde = A0/(Omega * sqrt(pi));
L       = 6.4e-3;        % 6400 um

% ---- 固定的延时（按你的要求）----
t_delay_x = 5088.85*6.4*1e-15;    % s

% ---------- Time arrays ----------
% 左侧波形用更致密的时间轴；右侧采样用于 Stokes
t3_time   = linspace(3.24/2*10^3*t0, 3.275/2*10^3*t0, 5000);
t3_values = linspace(3.250/2*10^3*t0, 3.270/2*10^3*t0, 101);
N         = numel(t3_values);

% ---------- Precompute common terms ----------
% Ey (左侧波形，用致密时间轴；右侧Stokes也复用相同物理表达式)
A_all   = 1/Omega^2 - 1i * k2 * L / 2;

% 左侧 Ey(t)
B_left   = k1 * L - t3_time;
Ey_left0 = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3_time) ...
           .* sqrt(pi ./ A_all) .* exp(-B_left.^2 ./ (4*A_all));
Ey_left  = 2 * real(Ey_left0);

% 右侧 Ey(t3_values)
B_vec     = k1 * L - t3_values;
Ey_right0 = E0tilde .* exp(1i*k*L) .* exp(-1i*omega0*t3_values) ...
            .* sqrt(pi./A_all) .* exp(-(B_vec.^2)./(4*A_all));

% ---------- Figure & layout ----------
fig = figure('Color','w','Position',[60 60 1350 600]);

% 1x9 网格：左 4，右 5 → 总体 4:5 布局
tiled = tiledlayout(fig,1,9,'Padding','compact','TileSpacing','compact');

% ===== Left panel: Ex/Ey (4/9 宽，纵横比 2:3) =====
ax1 = nexttile(tiled,[1 4]); hold(ax1,'on'); grid(ax1,'on');
ax = gca;
ax.Color = 'w';                % 背景白
ax.GridColor = [0.6 0.6 0.6];  % 深灰网格
ax.GridAlpha = 0.6;            % 半透明
ax.LineWidth = 0.8;              % 加粗坐标轴
set(ax1,'FontSize',14);
xlabel(ax1,'time [fs]','Interpreter','latex','FontSize',18);
ylabel(ax1,'$E(t)$ [a.u.]','Interpreter','latex','FontSize',18);
pbaspect(ax1,[3 2 1]);  % 宽:高 = 3:2

% Ex(t)（包含高斯包络 + 载频相位；不乘 e^{ik0L}，只影响全局相位）
Ex_left_c = A0 * exp(-(t3_time - t_delay_x).^2 / t0^2) .* exp(-1i*omega0*t3_time);
Ex_left   = 2 * real(Ex_left_c);  % 取实部并放大到与 Ey 同尺度可视


hEy = plot(ax1, t3_time*1e15, Ey_left, 'Color',[0   0   1],   'LineWidth',1.2);
hEx = plot(ax1, t3_time*1e15, Ex_left, 'Color',[0.5 0.5 0.5], 'LineWidth',1.2);
%legend(ax1, {'$E_y$','$E_x$'}, 'Interpreter','latex','FontSize',16,'Location','best');
legend(ax1, [hEx, hEy], {'$E_x$', '$E_y$'}, ...
       'Interpreter','latex','FontSize',16,'Location','best');

% ===== Right panel: Poincaré sphere (5/9 宽) =====
ax2 = nexttile(tiled,[1 5]); hold(ax2,'on');
view(ax2,120,30);
axis(ax2,'equal'); axis(ax2,[-1 1 -1 1 -1 1]); grid(ax2,'on');
set(ax2,'FontSize',14);
xticks(ax2,[-1 0 1]); yticks(ax2,[-1 0 1]); zticks(ax2,[-1 0 1]);
xlabel(ax2,'$s_1$','FontSize',26,'Interpreter','latex');
ylabel(ax2,'$s_2$','FontSize',26,'Interpreter','latex');
zlabel(ax2,'$s_3$','FontSize',26,'Interpreter','latex');

ax = gca;
grid on
ax.Color = 'w';                % 背景白
ax.GridColor = [0.6 0.6 0.6];  % 深灰网格
ax.GridAlpha = 0.6;            % 半透明
ax.LineWidth = 1;              % 加粗坐标轴

% 彩色四象限球
[xs, ys, zs] = sphere(100);
colors4 = [1   0.6 0.6;   % Q1 (y>=0, z>=0)
           0.6 0.6 1;     % Q2 (y<0,  z>=0)
           0.6 1   0.6;   % Q3 (y<0,  z<0)
           1   1   0.6];  % Q4 (y>=0, z<0)
C = zeros([size(ys),3]);
for i = 1:size(ys,1)
    for j = 1:size(ys,2)
        y = ys(i,j); z = zs(i,j);
        if y >= 0 && z >= 0
            C(i,j,:) = colors4(1,:);
        elseif y < 0 && z >= 0
            C(i,j,:) = colors4(2,:);
        elseif y < 0 && z < 0
            C(i,j,:) = colors4(3,:);
        else
            C(i,j,:) = colors4(4,:);
        end
    end
end
surf(ax2, xs, ys, zs, C, 'EdgeColor','none', 'FaceAlpha',0.15, 'HandleVisibility','off');
plot3(ax2,[-1 1],[0 0],[0 0],'--','Color',[0.4 0.4 0.4],'LineWidth',0.8,'HandleVisibility','off');
plot3(ax2,[0 0],[-1 1],[0 0],'--','Color',[0.4 0.4 0.4],'LineWidth',0.8,'HandleVisibility','off');
plot3(ax2,[0 0],[0 0],[-1 1],'--','Color',[0.4 0.4 0.4],'LineWidth',0.8,'HandleVisibility','off');
text(ax2,  1.1, 0,  -0.1, 'I', 'FontSize',20,'Interpreter','latex','HorizontalAlignment','center');
text(ax2, -1.1, 0.1, 0,    'F', 'FontSize',20,'Interpreter','latex','HorizontalAlignment','center');

% ---- Stokes on sphere (使用上述固定 t_delay_x) ----
Ex_r = A0 .* exp(-((t3_values - t_delay_x).^2)/t0^2) .* exp(-1i*omega0*t3_values) .* exp(1i*k0*L);
Ey_r = Ey_right0;

Ex2 = abs(Ex_r).^2; 
Ey2 = abs(Ey_r).^2;
S0  = max(Ex2 + Ey2, eps);
S1  = Ex2 - Ey2;
S2  = 2*real(Ex_r .* conj(Ey_r));
S3  = -2*imag(Ex_r .* conj(Ey_r));
s1n = S1 ./ S0; s2n = S2 ./ S0; s3n = S3 ./ S0;

% 样条平滑轨迹 + 彩色散点
t_idx = 1:numel(t3_values);
tt    = linspace(1, numel(t3_values), 20*numel(t3_values));
s1_spline = spline(t_idx, s1n, tt);
s2_spline = spline(t_idx, s2n, tt);
s3_spline = spline(t_idx, s3n, tt);

plot3(ax2, s1_spline, s2_spline, s3_spline, 'Color',[0 0 1], 'LineWidth', 2);
scatter3(ax2, s1n, s2n, s3n, 35, linspace(0,1,numel(t3_values)), ...
    'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.85,'HandleVisibility','off');
colormap(ax2,'jet'); clim(ax2,[0 1]);

% ========= 总标题 =========
%title(tiled, sprintf('$\\Delta t = %.2f\\,\\mathrm{fs}$', t_delay_x*1e15),'Interpreter','latex','FontSize',24);

% ---------- Export ----------
exportgraphics(gcf, 'Poincare_with_pulse.png', 'Resolution', 300);
disp('Saved figure: Poincare_with_pulse.png');

%% 17 200312 Ex Ey Pulse waveform only (Ex + Ey)

% ---------- Constants ----------
c = 3e8;
lambda0 = 800e-9;
omega0  = 2*pi*c / lambda0;
k0      = 2*pi / lambda0;
n       = 1.51078;

% ---------- Dispersion ----------
dn_dlambda   = -0.0198418e6;       % /m
d2n_dlambda2 = 0.0492482e12;       % /m^2
v_phi = c / n;
k     = omega0 / v_phi;
k1    = (1/c) * (n - lambda0 * dn_dlambda);
k2    = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

% ---------- Pulse / scan params ----------
t0      = 2e-14;         % 20 fs
Omega   = 2 / t0;
A0      = 1/2;
E0tilde = A0/(Omega * sqrt(pi));
L       = 6.4e-3;        % 6400 um

% ---- 固定延时 ----
t_delay_x = 5088.85*6.4*1e-15;

% ---------- Time arrays ----------
t3_time = linspace(3.24/2*10^3*t0, 3.275/2*10^3*t0, 5000);

% ---------- Precompute ----------
A_all   = 1/Omega^2 - 1i * k2 * L / 2;

B_left   = k1 * L - t3_time;
Ey_left0 = E0tilde .* exp(1i * k * L) .* exp(-1i * omega0 * t3_time) ...
           .* sqrt(pi ./ A_all) .* exp(-B_left.^2 ./ (4*A_all));
Ey_left  = 2 * real(Ey_left0);

% Ex pulse
Ex_left_c = A0 * exp(-(t3_time - t_delay_x).^2 / t0^2) .* exp(-1i*omega0*t3_time);
Ex_left   = 2 * real(Ex_left_c);

% ---------- Plot ----------
figure('Color','w','Position',[200 200 800 500])
hold on
grid on

ax = gca;
ax.Color = 'w';
ax.GridColor = [0.6 0.6 0.6];
ax.GridAlpha = 0.6;
ax.LineWidth = 1;

plot(t3_time*1e15, Ex_left,'Color',[0.5 0.5 0.5],'LineWidth',1.2)
plot(t3_time*1e15, Ey_left,'b','LineWidth',1.2)

xlabel('time [fs]','Interpreter','latex','FontSize',18)
ylabel('$E(t)$ [a.u.]','Interpreter','latex','FontSize',18)

legend({'$E_x$','$E_y$'},'Interpreter','latex','FontSize',16,'Location','best')

set(gca,'FontSize',14)

% ---------- Export ----------
exportgraphics(gcf,'pulse_waveform.png','Resolution',300);
disp('Saved figure: pulse_waveform.png');

%% 18  Pulse intensity envelope I = |Ex|^2 + |Ey|^2

% ---------- Constants ----------
c = 3e8;
lambda0 = 800e-9;
omega0  = 2*pi*c / lambda0;
k0      = 2*pi / lambda0;
n       = 1.51078;

% ---------- Dispersion ----------
dn_dlambda   = -0.0198418e6;       
d2n_dlambda2 = 0.0492482e12;       

v_phi = c / n;
k     = omega0 / v_phi;

k1 = (1/c) * (n - lambda0 * dn_dlambda);
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

% ---------- Pulse params ----------
t0      = 2e-14;      % 20 fs
Omega   = 2 / t0;
A0      = 1/2;
E0tilde = A0/(Omega * sqrt(pi));
L       = 6.4e-3;

t_delay_x = 5088.85*6.4*1e-15;

% ---------- Time ----------
t3_time = linspace(3.24/2*10^3*t0, 3.275/2*10^3*t0, 5000);

% ---------- Ey (complex field) ----------
A_all   = 1/Omega^2 - 1i * k2 * L / 2;
B_left  = k1 * L - t3_time;

Ey_left0 = E0tilde .* exp(1i*k*L) .* exp(-1i*omega0*t3_time) ...
           .* sqrt(pi./A_all) .* exp(-B_left.^2./(4*A_all));

% ---------- Ex (complex field) ----------
Ex_c = A0 * exp(-(t3_time - t_delay_x).^2 / t0^2) ...
      .* exp(-1i*omega0*t3_time);

% ---------- Intensity ----------
I = abs(Ex_c).^2 + abs(Ey_left0).^2;

% ---------- Plot ----------
figure('Color','w','Position',[200 200 800 500])
hold on
grid on

ax = gca;
ax.Color = 'w';
ax.GridColor = [0.6 0.6 0.6];
ax.GridAlpha = 0.6;
ax.LineWidth = 1;

plot(t3_time*1e15, I,'k','LineWidth',1.8)

xlabel('time [fs]','Interpreter','latex','FontSize',18)
ylabel('$I(t)$ [a.u.]','Interpreter','latex','FontSize',18)

set(gca,'FontSize',14)

% ---------- Export ----------
exportgraphics(gcf,'pulse_intensity_envelope.png','Resolution',300);
disp('Saved figure: pulse_intensity_envelope.png');

%% 19 Pulse intensity envelope for three time delays

% ---------- Constants ----------
c = 3e8;
lambda0 = 800e-9;
omega0  = 2*pi*c / lambda0;
k0      = 2*pi / lambda0;
n       = 1.51078;

% ---------- Dispersion ----------
dn_dlambda   = -0.0198418e6;
d2n_dlambda2 = 0.0492482e12;

v_phi = c / n;
k     = omega0 / v_phi;

k1 = (1/c) * (n - lambda0 * dn_dlambda);
k2 = (lambda0^3 / (2*pi*c^2)) * d2n_dlambda2;

% ---------- Pulse params ----------
t0      = 2e-14;
Omega   = 2 / t0;
A0      = 1/2;
E0tilde = A0/(Omega * sqrt(pi));
L       = 6.4e-3;

% ---------- Time delays (fs -> s) ----------
delays_fs = [32548 32568 32588];
delays    = delays_fs * 1e-15;

% ---------- Time axis ----------
t3_time = linspace(3.24/2*10^3*t0, 3.275/2*10^3*t0, 5000);

% ---------- Ey (complex field) ----------
A_all   = 1/Omega^2 - 1i * k2 * L / 2;
B_left  = k1 * L - t3_time;

Ey_left0 = E0tilde .* exp(1i*k*L) .* exp(-1i*omega0*t3_time) ...
           .* sqrt(pi./A_all) .* exp(-B_left.^2./(4*A_all));

% ---------- Plot ----------
figure('Color','w','Position',[200 200 850 520])
hold on
grid on

ax = gca;
ax.Color = 'w';
ax.GridColor = [0.6 0.6 0.6];
ax.GridAlpha = 0.6;
ax.LineWidth = 1;

colors = lines(3);

for i = 1:3
    
    % Ex for each delay
    Ex_c = A0 * exp(-(t3_time - delays(i)).^2 / t0^2) ...
          .* exp(-1i*omega0*t3_time);

    % Intensity
    I = abs(Ex_c).^2 + abs(Ey_left0).^2;

    % Plot
    plot(t3_time*1e15, I, 'LineWidth',2,'Color',colors(i,:))

end

xlabel('time [fs]','Interpreter','latex','FontSize',18)
ylabel('$I(t)$ [a.u.]','Interpreter','latex','FontSize',18)

legend({'32548 fs','32568 fs','32588 fs'}, ...
       'Location','best','FontSize',14)

set(gca,'FontSize',14)

% ---------- Export ----------
exportgraphics(gcf,'pulse_intensity_three_delays.png','Resolution',300);
disp('Saved figure: pulse_intensity_three_delays.png');


