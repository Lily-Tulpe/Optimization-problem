% 参数初始化
Pc = 1; % Ps1 is the sensing power, Pc is the communication power, Ps1 = Pc < 1W
hc = 0.5; % hs1^2 is the scalar reflection coefficient, hc is the channel coefficient, hs1 < hc < 1
d = 480; % d is the length of packet, d = 2^n * 10^m
sigma2 = 0.01; % sigma^2 is the noise power, sigma^2 = 0.01
Ts = 0.02; % Ts is the duration of a single symbol, Ts < 1ms, Ts = 0.025ms / 0.02ms
Tmax = 5; % Tmax is the total time, Tmax < 100ms, Tmax < 10ms, Tmax < 5ms
 
gamma = (Pc*(hc^2))/(sigma2); 
V = 1-1/((1+gamma)^2); 
C = log2(1+gamma); 

Ps1 = 1.2;
hs1 = 0.15;
delta1 = 0.0001; % delta1 is the probability of false alarm (PFA) threshold, 1e-5 < delta1 < 1e-3

Ps2 = 1;
hs2 = 0.1;
delta2 = 0.0001;

err_s_th = 0.01;
err_c_th = 0.01;

% 感知错误率阈值
y = qfuncinv(err_s_th);
Ls_sqrt = (y*sqrt(2*Ps1*sigma2*(hs1^2))+sqrt((y.^2)*2*Ps1*sigma2*(hs1^2)+4*Ps1*(hs1^2)*sigma2*(-log(delta1))))/(2*Ps1*(hs1^2));
Ls = ceil(Ls_sqrt.^2);
ts1_th = Ls*Ts;

Ls_sqrt = (y*sqrt(2*Ps2*sigma2*(hs2^2))+sqrt((y.^2)*2*Ps2*sigma2*(hs2^2)+4*Ps2*(hs2^2)*sigma2*(-log(delta2))))/(2*Ps2*(hs2^2));
Ls = ceil(Ls_sqrt.^2);
ts2_th = Ls*Ts;

% 通信错误率阈值
y = qfuncinv(err_c_th);
Lc_sqrt = (y/log(2)+sqrt((y/log(2)).^2+4*C*d/V))/(2*C/sqrt(V));
Lc = ceil(Lc_sqrt.^2);
tc_th = Lc*Ts;

% 黄金分割法参数
tau = (sqrt(5) - 1) / 2; % 黄金分割比例
tol = 1e-6; % 收敛精度
max_iter = 100; % 最大迭代次数

a = [ts1_th; Tmax-tc_th-ts1_th]; % 下界
b = [Tmax-tc_th-ts2_th; ts2_th]; % 上界

% 黄金分割法
for iter = 1:max_iter

    x1 = a + (1 - tau) * (b - a);
    x2 = a + tau * (b - a);
    
    % 计算 err_s
    err_s_x1 = calculate_error_s(x1, Ps1, Ps2, hs1, hs2, sigma2, delta1, delta2, Ts);
    
    err_s_x2 = calculate_error_s(x2, Ps1, Ps2, hs1, hs2, sigma2, delta1, delta2, Ts);
                   
    % 更新区间
    if err_s_x1 < err_s_x2
        b = x2;            
    else
        a = x1;
    end
    
    % 检查收敛
    if norm(b - a) < tol
        break;
    end
end    

% 输出结果
x_opt = (a+b)/2;
err_s_opt = calculate_error_s(x_opt, Ps1, Ps2, hs1, hs2, sigma2, delta1, delta2, Ts);

fprintf('Optimal tc=%.4f, ts1=%.4f, ts2=%.4f\n', tc_th, x_opt(1), x_opt(2));
fprintf('Optimal err_s=%e\n',err_s_opt);

% 计算感知错误率的函数
function err_s = calculate_error_s(x, Ps1, Ps2, hs1, hs2, sigma2, delta1, delta2, Ts)
    Ls1 = x(1) / Ts;
    err_s1 = qfunc((Ps1 * Ls1 * hs1^2 - (sigma2) * (-log(delta1))) ./ (sqrt(2 * Ps1 * Ls1 * (sigma2) * hs1^2)));
    
    Ls2 = x(2) / Ts;
    err_s2 = qfunc((Ps2 * Ls2 * hs2^2 - (sigma2) * (-log(delta2))) ./ (sqrt(2 * Ps2 * Ls2 * (sigma2) * hs2^2)));
     
    err_s = err_s1 + err_s2 - err_s1 .* err_s2;
end