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
title(['Electric Field Distribution in Frequency Domain (T = ' num2str(T) ' s)']);
grid on;

%%
% 清理环境
clear; clc; close all;

% 设定参数
E0 = 1;                        % 振幅
w0 = 2.9767e15;                % 中心角频率 (Hz)
T = 15e-15;                     % 时间宽度 (s)
Omega = 2/T;                   % 频率带宽 (Hz)
sigma = Omega / sqrt(2);       % 高斯分布的标准差

% 定义频率范围
w = linspace(w0 - 3*Omega, w0 + 3*Omega, 1000); 

% 计算频谱中的电场分布（修正后的高斯分布）
E_w = E0 * exp(-((w - w0).^2) / (2*sigma^2));

% 计算时间域电场（使用傅里叶变换）
t = linspace(-5*T, 5*T, 1000); % 时间范围
E_t = fftshift(ifft(ifftshift(E_w))); % 逆傅里叶变换并对齐中心

% 画出频域电场
figure;

plot(w, abs(E_w), 'b', 'LineWidth', 2);
xlabel('Frequency \omega (Hz)');
ylabel('|{E}(\omega)|');
title('Electric Field in Frequency Domain');
grid on;
xlim([w0 - 3*Omega, w0 + 3*Omega]); % 保证w0为中心

%%
% 设定参数
E0 = 1;                        % 振幅
w0 = 2.9767e15;                % 中心角频率 (Hz)
T_values = [1e-15, 2e-15, 3e-15, 10e-15];  % 不同的时间宽度值

colors = ['r', 'c','b', 'g'];  % 颜色定义
figure;

% 频率范围
w_min = w0 - 3 * max(2 ./ T_values);
w_max = w0 + 3 * max(2 ./ T_values);
w = linspace(w_min, w_max, 1000);

% 频域电场

hold on;
for i = 1:length(T_values)
    T = T_values(i);
    Omega = 2 / T;
   
    % 计算频谱中的电场分布（修正后的高斯分布）
    E_w = E0 * exp(-((w - w0).^2) / (1 * Omega^2));

    plot(w, abs(E_w), 'Color', colors(i), 'LineWidth', 2, 'DisplayName', ['T = ', num2str(T*1e15), ' fs']);
end

xlabel('Frequency \omega (Hz)');
ylabel('|{E}(\omega)|');
title('Electric Field in Frequency Domain');
legend show;
grid on;
xlim([w_min, w_max]);
hold off;

% 时间范围
t = linspace(-5 * max(T_values), 5 * max(T_values), 1000);
