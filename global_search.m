% 参数初始化
Pc = 1; % Ps is the sensing power,Pc is the communication power,Ps=Pc<1W
hc = 0.5; % hs^2 is the scalar refection coefficient,hc is the channel coefficient,hs<hc<1
sigma2 = 0.01; % sigma^2 is the noise power,sigma^2=0.01
d = 480; % d is the length of packet,d=2^n*10^m
Ts = 0.02; % Ts is the duration of a single symbol,Ts<1ms,Ts=0.025ms/0.02ms
Tmax = 5; % Tmax is the total time,Tmax<100ms,Tmax<10ms,Tmax<5ms

hs1 = 0.15;
Ps1 = 1.2;
delta1 = 0.0001; % delta is the probability of false alarm (PFA) threshold,1e-5<delta<1e-3

hs2 = 0.1;
Ps2 = 1;
delta2 = 0.0001;

err_s_th = 0.01;
err_c_th = 0.01;

% 感知时间 ts1 相关计算
ts1 = [Ts:Ts:Tmax];
Ls1 = ts1/Ts;
err_s1 = qfunc((Ps1*Ls1*hs1^2-(sigma2)*(-log(delta1)))./(sqrt(2*Ps1*Ls1*(sigma2)*hs1^2))); 
ts1(err_s1 > err_s_th) = NaN;

% 感知时间 ts2 相关计算
ts2 = [Ts:Ts:Tmax];
Ls2 = ts2/Ts;
err_s2 = qfunc((Ps2*Ls2*hs2^2-(sigma2)*(-log(delta2)))./(sqrt(2*Ps2*Ls2*(sigma2)*hs2^2))); 
ts2(err_s2 > err_s_th) = NaN;

% 通信时间 tc 相关计算
gamma = (Pc*(hc^2))/(sigma2); 
V = 1 - 1/((1 + gamma)^2); 
C = log2(1 + gamma); 

tc = [Ts:Ts:Tmax];
Lc = tc/Ts;
err_c = qfunc(sqrt(Lc/V).*(C - d./Lc)*log(2));
tc(err_c > err_c_th) = NaN;

% 生成三维网格点
[TC, TS1, TS2] = meshgrid(tc, ts1, ts2);
         
% 筛选出满足 tc + ts1 + ts2 ≤ Tmax 的点
valid_indices = ~isnan(TC) & ~isnan(TS1) & ~isnan(TS2) & (TC + TS1 + TS2 <= Tmax);

% 提取可行域的点
valid_tc = TC(valid_indices);
valid_ts1 = TS1(valid_indices);
valid_ts2 = TS2(valid_indices);

% 计算可行域内目标函数值
Ls1 = valid_ts1 / Ts;
err_s1 = qfunc((Ps1 * Ls1 * hs1^2 - (sigma2) * (-log(delta1))) ./ (sqrt(2 * Ps1 * Ls1 * (sigma2) * hs1^2)));

Ls2 = valid_ts2 / Ts;
err_s2 = qfunc((Ps2 * Ls2 * hs2^2 - (sigma2) * (-log(delta2))) ./ (sqrt(2 * Ps2 * Ls2 * (sigma2) * hs2^2)));

err_s = err_s1 + err_s2 - err_s1 .* err_s2;

% 求出最优解
[err_s_opt,I]=min(err_s);

tc_opt = valid_tc(I);
ts1_opt = valid_ts1(I);
ts2_opt = valid_ts2(I);

fprintf('最优点为：tc=%.4f, ts1=%.4f, ts2=%.4f\n', tc_opt, ts1_opt, ts2_opt);
fprintf('最优解为：err_s=%e\n',err_s_opt);

