clear; clc;
%% 前期定义
% 定义各参数大小
R0 = 5e-9; %光斑半径
U0 = 0.5; %电子初始动能
U = 12.67e3; %施加电压大小
alpha = 10.2; %电子最大出射角度
n = 10000; %计算步数
L = 1e-3; %运动距离
J = 7e-4; %电流密度
e0 = 8.854e-12; %真空介电常数
e = 1.60217662e-19;
m = 9.1093837e-31;

% 衍生函数大小
B = B_calculation(U, L, U0, alpha, 1); %需要施加的磁场大小
% T = 2*pi*m/(e*B);
T =L*m/(e*U)*sqrt(2*e*U0/m)*(sqrt(1+U/(U0*cosd(alpha)^2))-1)*cosd(alpha); %总运动时间，以及粒子旋转的周期
dt = T/n; %时间步长

%% 离散求和方式求影响
% 预计算公共变量
i_vec = 0:n;              % 生成索引向量 [0,1,...,n]
ti = i_vec * dt;          % 时间向量

% 向量化计算 R, t, h, E (消除第一个循环)
R = R0 + 2/B*sqrt(2*m*U0/e)*sind(alpha)*sin(e*B/(2*m)*ti);
t = ti;
h = sqrt(2*e*U0/m)*sind(alpha).*ti + 0.5*U*e/(L*m)*ti.^2;
E = J*pi*R0^2 ./ (2*pi*e0*R .* (sqrt(2*e*U0/m)*cosd(alpha) + U*e/(L*m)*ti));

% 向量化计算 dd1
cumsum_E = cumsum(E);     % 预计算累积和
sum_E_part(1:2) = 0;
sum_E_part(3:n+1) = cumsum_E(1:n-1);  % 对n≥3，sum_{k=1}^{n-2} = cumsum_E(n-2)
E_part = 0;
E_part(2:n+1) = E(1:n);

dd1 = (e/m) * (sum_E_part + 0.5*E_part) * dt^2;
d1 = sum(dd1);

%% 积分方式求影响
syms t1;
E1 = J*pi*R0^2 ./ (2*pi*e0*(R0 + 2/B*sqrt(2*m*U0/e)*sind(alpha)*sin(e*B/(2*m)*t1)) .* (sqrt(2*e*U0/m)*cosd(alpha) + U*e/(L*m)*t1));  % 示例函数

% 修正积分表达式：使用数值T_final作为上限
A_integral = (e/m) * int((T - t1)*E1, t1, 0, T);

% 离散求和部分修正：显式转换符号表达式为数值
E_discrete = double(subs(E1, t1, t));  % 强制转换为数值数组

B_sum = 0.5 * (e/m) * dt^2 * sum(E_discrete);

d2 = double(A_integral) + B_sum;  % 分步转换确保类型匹配

%% 结果展示
% 包络形状绘制
figure(2)
plot(h,R)

% 求和与积分结果展示
disp(['总位移 求和计算结果d1 = ', num2str(d1), '，', '积分计算结果d2 = ', num2str(d2) ]);

%% 绘制横向电场随位置变化
figure(3)
plot(h, E, 'b-', 'LineWidth', 2)
xlabel('轴向位置 h (m)')
ylabel('横向电场 E (V/m)')
title('电子束包络上的横向电场分布')
grid on

% 可选：添加对数坐标来更好地观察变化趋势
figure(4)
subplot(2,1,1)
semilogy(h, abs(E), 'r-', 'LineWidth', 2)
xlabel('轴向位置 h (m)')
ylabel('|E| (V/m)')
title('横向电场绝对值（对数坐标）')
grid on

subplot(2,1,2)
plot(h, E, 'b-', 'LineWidth', 2)
xlabel('轴向位置 h (m)')
ylabel('E (V/m)')
title('横向电场（线性坐标）')
grid on

% 显示电场统计信息
fprintf('横向电场统计信息:\n');
fprintf('最大值: %.4e V/m\n', max(E));
fprintf('最小值: %.4e V/m\n', min(E));
fprintf('平均值: %.4e V/m\n', mean(E));
fprintf('电场变化范围: %.4e V/m\n', range(E));

% 详细电场分析
figure(5)
subplot(2,2,1)
plot(h, E, 'b-', 'LineWidth', 1.5)
xlabel('轴向位置 h (m)')
ylabel('E (V/m)')
title('横向电场分布')
grid on

subplot(2,2,2)
histogram(E, 50, 'FaceColor', 'green', 'FaceAlpha', 0.7)
xlabel('电场强度 (V/m)')
ylabel('频数')
title('电场强度分布直方图')

subplot(2,2,3)
plot(t, E, 'm-', 'LineWidth', 1.5)
xlabel('时间 t (s)')
ylabel('E (V/m)')
title('横向电场随时间变化')
grid on

subplot(2,2,4)
plot(h, R, 'k-', 'LineWidth', 1.5)
hold on
plot(h, abs(E)/max(abs(E))*max(R), 'r--', 'LineWidth', 1.5)
xlabel('轴向位置 h (m)')
ylabel('半径/归一化电场')
title('束包络与归一化电场对比')
legend('束半径 R', '归一化|E|')
grid on