%% 定义参数
A0 = 1; % 振幅
T = 2e-15; % 设定时间尺度参数，假设T为1 fs
omega0 = 2.9767e15; % 角频率 (rad/s)

% 定义时间范围
t = linspace(-5*T, 5*T, 1000); % 在较宽的范围内取点

% 计算复高斯脉冲
E_t = A0 * exp(-t.^2 / T^2) .* exp(-1i * omega0 * t);

% 计算实部和包络
E_real = real(E_t);
Envelope = abs(A0 * exp(-t.^2 / T^2));

% 绘制图像
figure;
plot(t, E_real, 'b', 'LineWidth', 1.5); hold on;
plot(t, Envelope, 'r--', 'LineWidth', 1.5); % 绘制包络
plot(t, -Envelope, 'r--', 'LineWidth', 1.5); % 绘制包络
grid on;
xlabel('Time (t) [s]');
ylabel('Amplitude');
title('Gaussian Pulse with Carrier Frequency');
legend('Re(E(0,t))', 'Envelope', 'Location', 'best');

%%  Gaussian Pulse with 375 THz Carrier (800 nm)
A0 = 1;                        % 振幅
T = 5e-15;                   % 脉冲宽度（100 fs）
f0 = 375e12;                   % 中心频率（375 THz）
omega0 = 2 * pi * f0;          % 角频率 (rad/s)

% 定义时间范围（5倍脉宽）
t = linspace(-5*T, 5*T, 1000); % 时间轴，覆盖足够脉冲范围

% 计算复高斯脉冲
E_t = A0 * exp(-t.^2 / T^2) .* exp(-1i * omega0 * t);

% 计算实部和包络
E_real = remiroal(E_t);
Envelope = abs(A0 * exp(-t.^2 / T^2));

% 绘图
figure;
plot(t, E_real, 'b', 'LineWidth', 1.5); hold on;
plot(t, Envelope, 'r--', 'LineWidth', 1.5);
plot(t, -Envelope, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Time (t) [s]');
ylabel('Amplitude');
title('Gaussian Pulse (5 fs) with 375 THz Carrier (800 nm)');
legend('Re(E(t))', 'Envelope', 'Location', 'best');


