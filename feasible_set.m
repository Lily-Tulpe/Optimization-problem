%%
% 定义符号变量
syms ts1 ts2 tc; 

% 参数初始化
hs = 0.1;
hc = 0.5; % hs^2 is the scalar refection coefficient,hr is the channel coefficient,hs<hr<1
sigma2 = 0.01; % sigma^2 is the noise power,sigma^2=0.01
Ps = 1;
Pc = 1; % Ps is the sensing power,Pr is the communication power,Ps=Pr<1W
d = 480; % d is the length of packet,d=2^n*10^m
Ts = 0.02; % Ts is the duration of a single symbol,Ts<1ms,Ts=0.025ms/0.02ms
Tmax = 5; % Tmax is the total time,Tmax<100ms,Tmax<10ms,Tmax<5ms
delta = 0.0001; % delta is the probability of false alarm (PFA) threshold,1e-5<delta<1e-3

err_s_th = 0.01;
err_c_th = 0.01;

n=5;
% 感知时间 ts1 相关计算
ts1 = [Ts:n*Ts:Tmax];
Ls1 = ts1/Ts;
err_s1 = qfunc((Ps*Ls1*hs^2-(sigma2)*(-log(delta)))./(sqrt(2*Ps*Ls1*(sigma2)*hs^2))); 
ts1(err_s1 >= err_s_th) = NaN;

% 感知时间 ts2 相关计算
ts2 = [Ts:n*Ts:Tmax];
Ls2 = ts2/Ts;
err_s2 = qfunc((Ps*Ls2*hs^2-(sigma2)*(-log(delta)))./(sqrt(2*Ps*Ls2*(sigma2)*hs^2))); 
ts2(err_s2 >= err_s_th) = NaN;

% 通信时间 tc 相关计算
gamma = (Pc*(hc^2))/(sigma2); 
V = 1 - 1/((1 + gamma)^2); 
C = log2(1 + gamma); 

tc = [Ts:n*Ts:Tmax];
Lc = tc/Ts;
err_c = qfunc(sqrt(Lc/V).*(C - d./Lc)*log(2));
tc(err_c >= err_c_th) = NaN;

% 生成三维网格点
[TC, TS1, TS2] = meshgrid(tc, ts1, ts2);

% 筛选出满足 tc + ts1 + ts2 ≤ Tmax 的点
valid_indices = ~isnan(TC) & ~isnan(TS1) & ~isnan(TS2) & (TC + TS1 + TS2 <= Tmax);

% 提取可行域的点
valid_tc_sampled = TC(valid_indices);
valid_ts1_sampled = TS1(valid_indices);
valid_ts2_sampled = TS2(valid_indices);

Ls1 = valid_ts1_sampled/Ts;
err_s1_sampled = qfunc((Ps*Ls1*hs^2-(sigma2)*(-log(delta)))./(sqrt(2*Ps*Ls1*(sigma2)*hs^2))); 
Ls2 = valid_ts2_sampled/Ts;
err_s2_sampled = qfunc((Ps*Ls2*hs^2-(sigma2)*(-log(delta)))./(sqrt(2*Ps*Ls2*(sigma2)*hs^2)));
err_s_sampled = err_s1_sampled+err_s2_sampled-err_s1_sampled.*err_s2_sampled;

% 绘制三维散点图，并根据目标函数的值染色
figure;
scatter3(valid_tc_sampled, valid_ts1_sampled, valid_ts2_sampled, 10, err_s_sampled, 'filled');
hold on;

%{
% 绘制平面 tc = 2.2
[ts1_plane, ts2_plane] = meshgrid(linspace(0, Tmax, 100), linspace(0, Tmax, 100));
tc_plane_tc = 2.2 * ones(size(ts1_plane));
surf(tc_plane_tc, ts1_plane, ts2_plane, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% 绘制平面 ts1 = 0.52
[tc_plane, ts2_plane] = meshgrid(linspace(0, Tmax, 100), linspace(0, Tmax, 100));
ts1_plane_ts1 = 0.52 * ones(size(tc_plane));
surf(tc_plane, ts1_plane_ts1, ts2_plane, 'FaceColor', 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% 绘制平面 ts2 = 0.52
[tc_plane, ts1_plane] = meshgrid(linspace(0, Tmax, 100), linspace(0, Tmax, 100));
ts2_plane_ts2 = 0.52 * ones(size(tc_plane));
surf(tc_plane, ts1_plane, ts2_plane_ts2, 'FaceColor', 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% 绘制平面 tc + ts1 + ts2 = Tmax
[tc_plane, ts1_plane] = meshgrid(linspace(0, Tmax, 100), linspace(0, Tmax, 100));
ts2_plane_total = Tmax - tc_plane - ts1_plane;
surf(tc_plane, ts1_plane, ts2_plane_total, 'FaceColor', 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
%}

% 添加颜色条
colorbar;
xlabel('tc');
ylabel('ts1');
zlabel('ts2');
ylim([0.5 2.3]);zlim([0.5 2.3]);
title('Feasible Region with Colormap (Sampled Points)');
hold off;
