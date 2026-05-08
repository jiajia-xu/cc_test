%% 1 Sellmeier Equation Plot for BK7 Glass

% Wavelength range in micrometers (μm)
lambda = linspace(0.2, 1.6, 500);  % 避免从0开始以防除以0

% Sellmeier 系数（BK7玻璃）
B1 = 1.03961212;
B2 = 0.231792344;
B3 = 1.01046945;
C1 = 6.00069867e-3;
C2 = 2.00179144e-2;
C3 = 1.03560653e2;

% 根据 Sellmeier 方程计算折射率
n_squared = 1 + ...
    (B1 * lambda.^2) ./ (lambda.^2 - C1) + ...
    (B2 * lambda.^2) ./ (lambda.^2 - C2) + ...
    (B3 * lambda.^2) ./ (lambda.^2 - C3);
n = sqrt(n_squared);

% 绘制折射率曲线
figure;
plot(lambda, n, 'r-', 'LineWidth', 2);
xlabel('Wavelength \lambda (\mum)', 'FontSize', 12);
ylabel('Refractive index n(\lambda)', 'FontSize', 12);
title('Sellmeier Equation for BK7 Glass', 'FontSize', 14);
grid on;
xlim([0.2 1.6]);
ylim([1.5 1.6]);

%% 2 Sellmeier Equation Plot for BK7 Glass

% Wavelength range in micrometers (μm)
lambda = linspace(0.2, 1.6, 500);  % 避免从0开始以防除以0

% Sellmeier 系数（BK7玻璃）
B1 = 1.03961212;
B2 = 0.231792344;
B3 = 1.01046945;
C1 = 6.00069867e-3;
C2 = 2.00179144e-2;
C3 = 1.03560653e2;

% 根据 Sellmeier 方程计算折射率
n_squared = 1 + ...
    (B1 * lambda.^2) ./ (lambda.^2 - C1) + ...
    (B2 * lambda.^2) ./ (lambda.^2 - C2) + ...
    (B3 * lambda.^2) ./ (lambda.^2 - C3);
n = sqrt(n_squared);

% 计算 0.8 μm 的折射率
lambda_point = 0.8;
n_point_sq = 1 + ...
    (B1 * lambda_point^2) / (lambda_point^2 - C1) + ...
    (B2 * lambda_point^2) / (lambda_point^2 - C2) + ...
    (B3 * lambda_point^2) / (lambda_point^2 - C3);
n_point = sqrt(n_point_sq);

% 绘制折射率曲线
figure;
plot(lambda, n, 'r-', 'LineWidth', 2);
hold on;

% 标记 0.8 μm 点
plot(lambda_point, n_point, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
text(lambda_point + 0.02, n_point, ...
    sprintf('(%0.1f \\mum, %.4f)', lambda_point, n_point), ...
    'FontSize', 12, 'Color', 'b');

xlabel('Wavelength \lambda (\mum)', 'FontSize', 12);
ylabel('Refractive index n(\lambda)', 'FontSize', 12);
title('Sellmeier Equation for BK7 Glass', 'FontSize', 14);
grid on;
xlim([0.2 1.6]);
ylim([1.5 1.6]);
%% 3
% Sellmeier Equation Plot for BK7 Glass

% Wavelength range in micrometers (μm)
lambda = linspace(0.2, 1.6, 500);  % 避免从0开始以防除以0

% Sellmeier 系数（BK7玻璃）
B1 = 1.03961212;
B2 = 0.231792344;
B3 = 1.01046945;
C1 = 6.00069867e-3;
C2 = 2.00179144e-2;
C3 = 1.03560653e2;

% 根据 Sellmeier 方程计算折射率
n_squared = 1 + ...
    (B1 * lambda.^2) ./ (lambda.^2 - C1) + ...
    (B2 * lambda.^2) ./ (lambda.^2 - C2) + ...
    (B3 * lambda.^2) ./ (lambda.^2 - C3);
n = sqrt(n_squared);

% 数值计算 dn/dλ
dn_dlambda = gradient(n, lambda);  % 使用梯度法进行数值微分

% 计算在 lambda = 0.8 μm 处的折射率和导数值
lambda_point = 0.8;
[~, idx] = min(abs(lambda - lambda_point));  % 找到最接近0.8的位置
n_point = n(idx);
dn_dlambda_point = dn_dlambda(idx);

% 绘图
figure;
plot(lambda, n, 'r-', 'LineWidth', 2);
hold on;

% 标记 0.8 μm 点
plot(lambda_point, n_point, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
text(lambda_point + 0.02, n_point, ...
    sprintf('(%0.1f \\mum, %.4f)', lambda_point, n_point), ...
    'FontSize', 12, 'Color', 'b');

% 显示 dn/dlambda 值
text(lambda_point + 0.02, n_point - 0.005, ...
    sprintf('dn/d\\lambda = %.4f /\\mum', dn_dlambda_point), ...
    'FontSize', 12, 'Color', 'k');

xlabel('Wavelength \lambda (\mum)', 'FontSize', 12);
ylabel('Refractive index n(\lambda)', 'FontSize', 12);
title('Sellmeier Equation for BK7 Glass', 'FontSize', 14);
grid on;
xlim([0.2 1.6]);
ylim([1.5 1.6]);

%% 4
% Sellmeier Equation Plot for BK7 Glass

% Wavelength range in micrometers (μm)
lambda = linspace(0.2, 1.6, 500);  % 避免从0开始以防除以0

% Sellmeier 系数（BK7玻璃）
B1 = 1.03961212;
B2 = 0.231792344;
B3 = 1.01046945;
C1 = 6.00069867e-3;
C2 = 2.00179144e-2;
C3 = 1.03560653e2;

% 根据 Sellmeier 方程计算折射率
n_squared = 1 + ...
    (B1 * lambda.^2) ./ (lambda.^2 - C1) + ...
    (B2 * lambda.^2) ./ (lambda.^2 - C2) + ...
    (B3 * lambda.^2) ./ (lambda.^2 - C3);
n = sqrt(n_squared);

% 数值计算一阶导数 dn/dλ 和二阶导数 d²n/dλ²
dn_dlambda = gradient(n, lambda);
d2n_dlambda2 = gradient(dn_dlambda, lambda);

% 计算在 lambda = 0.8 μm 处的折射率和导数值
lambda_point = 0.8;
[~, idx] = min(abs(lambda - lambda_point));  % 找到最接近0.8的位置
n_point = n(idx);
dn_dlambda_point = dn_dlambda(idx);
d2n_dlambda2_point = d2n_dlambda2(idx);

% 绘图
figure;
plot(lambda, n, 'r-', 'LineWidth', 2);
hold on;

% 标记 0.8 μm 点
plot(lambda_point, n_point, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
text(lambda_point + 0.02, n_point, ...
    sprintf('(%0.1f \\mum, %.4f)', lambda_point, n_point), ...
    'FontSize', 12, 'Color', 'b');

% 显示一阶导数值
text(lambda_point + 0.02, n_point - 0.005, ...
    sprintf('dn/d\\lambda = %.4f /\\mum', dn_dlambda_point), ...
    'FontSize', 12, 'Color', 'k');

% 显示二阶导数值
text(lambda_point + 0.02, n_point - 0.010, ...
    sprintf('d^2n/d\\lambda^2 = %.4f /\\mum^2', d2n_dlambda2_point), ...
    'FontSize', 12, 'Color', 'm');

xlabel('Wavelength \lambda (\mum)', 'FontSize', 12);
ylabel('Refractive index n(\lambda)', 'FontSize', 12);
title('Sellmeier Equation for BK7 Glass', 'FontSize', 14);
grid on;
xlim([0.2 1.6]);
ylim([1.5 1.6]);

