%%

% 常量
c = 3e8; % 光速 (m/s)
lambda0 = 0.8e-6; % 中心波长 (m)
n = 1.5108;

% 导数
dn_dlambda = -0.0198e6; % 转换为 /m
d2n_dlambda2 = 0.0491e12; % 转换为 /m^2

% ω
omega = 2 * pi * c / lambda0;

% 使用第二个公式：
term = (-lambda0^2 / (2 * pi * c)) * dn_dlambda;
dk_domega = (1 / c) * (omega * term + n);

% 输出结果
fprintf('dk/dω = %.6e s/m (λ₀)\n', dk_domega);


%%
% 常量
c = 3e8; % 光速 (m/s)
lambda0 = 0.8e-6; % 中心波长 (m)
n = 1.5108;
dn_dlambda = -0.0198e6; % 转换为 /m

% ω
omega = 2 * pi * c / lambda0;

% 计算分子和分母
numerator = (-2 * pi * c) / (n * omega) * dn_dlambda;
denominator = 1 + (2 * pi * c) / (omega * n^2) * dn_dlambda;

% 总表达式
dk_domega_1 = (1 / c) * (numerator / denominator + n);

% 输出结果
fprintf('dk/dω = %.6e s/m\n', dk_domega_1);
