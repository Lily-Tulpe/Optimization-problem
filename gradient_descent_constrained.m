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

% 初始化参数
max_iter = 100000; % 最大迭代次数
tol = 1e-6; % 收敛容差
c = 0.1; % Armijo 条件常数
rho = 0.5; % 步长缩减因子
alpha = 1; % 初始步长

% 初始值
t = [tc_th; ts1_th; ts2_th];
     
% 梯度下降迭代
for iter = 1:max_iter

    % 计算梯度
    grad_tc = 0;

    Ls1=t(2)/Ts;
    ws1=(Ps1*Ls1*hs1^2-(sigma2)*(-log(delta1)))./(sqrt(2*Ps1*Ls1*(sigma2)*hs1^2));
    err_s1 = qfunc(ws1);
    
    Ls2 = t(3)/Ts;
    ws2=(Ps2*Ls2*hs2^2-(sigma2)*(-log(delta2)))./(sqrt(2*Ps2*Ls2*(sigma2)*hs2^2));
    err_s2 = qfunc(ws2); 

    grad_ts1 = (-1/(2*sqrt(2*pi)*Ts))*exp(-(ws1.^2)/2).*((Ps1*hs1^2*sqrt(Ls1)+sigma2*(-log(delta1))*(1./sqrt(Ls1)))./(sqrt(2*Ps1*sigma2*hs1^2)*Ls1)).*(1-err_s2);   
    grad_ts2 = (-1/(2*sqrt(2*pi)*Ts))*exp(-(ws2.^2)/2).*((Ps2*hs2^2*sqrt(Ls2)+sigma2*(-log(delta2))*(1./sqrt(Ls2)))./(sqrt(2*Ps2*sigma2*hs2^2)*Ls2)).*(1-err_s1);
    
    gradient = [grad_tc; grad_ts1; grad_ts2];
 
    % 线搜索
    while true
        % 更新变量
        t_new = t - alpha * gradient;
        
        % 投影到约束空间
        t_new(1) = max(t_new(1), tc_th); % tc >= tc_th
        t_new(2) = max(t_new(2), ts1_th); % ts1 >= ts1_th
        t_new(3) = max(t_new(3), ts2_th); % ts2 >= ts2_th
            
        % 检查约束条件 tc + ts1 + ts2 <= Tmax
        if sum(t_new) > Tmax
            % 投影到约束条件
            t_new = Tmax*(t_new / sum(t_new));
        end
        
        % 计算目标函数值
        err_s = calculate_error_s(t, Ps1, hs1, Ps2, hs2, sigma2, delta1, delta2, Ts);
        err_s_new = calculate_error_s(t_new, Ps1, hs1, Ps2, hs2, sigma2, delta1, delta2, Ts);
        
        % 检查 Armijo 条件
        if err_s_new <= err_s - c * alpha * norm(gradient)^2
            break;
        end
        
        % 缩减步长
        alpha = rho * alpha;
    end
 
    % 检查收敛
    if norm(t_new-t) < tol
        break;
    end
    
    % 更新变量
    t = t_new;
end

err_s_opt = calculate_error_s(t, Ps1, hs1, Ps2, hs2, sigma2, delta1, delta2, Ts);
% 输出结果
fprintf('最优解: tc = %.4f, ts1 = %.4f, ts2 = %.4f\n', t(1), t(2), t(3));
fprintf('目标函数值: %e\n', err_s_opt);

function err_s = calculate_error_s(x, Ps1, hs1, Ps2, hs2, sigma2, delta1, delta2, Ts)
    Ls1 = x(2) / Ts;
    err_s1 = qfunc((Ps1 * Ls1 * hs1^2 - sigma2 * (-log(delta1))) ./ (sqrt(2 * Ps1 * Ls1 * sigma2 * hs1^2)));
    
    Ls2 = x(3) / Ts;
    err_s2 = qfunc((Ps2 * Ls2 * hs2^2 - sigma2 * (-log(delta2))) ./ (sqrt(2 * Ps2 * Ls2 * sigma2 * hs2^2)));
    
    err_s = err_s1 + err_s2 - err_s1 .* err_s2;
end
