c0 = 0.5;
T = 10.0;
tau = 0.2;

[t, Cn, En] = simulate_stage3(c0, T, tau);
C_exact = exact_sol(t, c0);
max_abs_error = max(abs(Cn - C_exact));
fprintf('tau = %.2f, maximum point-wise absolute error = %.4e\n',tau,max_abs_error);

figure('Position',[100,100,600,400]);
plot(t,Cn,'DisplayName','RK4 numerical');
hold on;
plot(t,C_exact,'--','DisplayName','Exact ODE solution');
xlabel('$t$','Interpreter','latex');
ylabel('$c(t)$','Interpreter','latex');
legend;
title(sprintf('RK4, $\\tau = %.2f$',tau),'Interpreter','latex');
hold off;

figure('Position',[100,100,600,400]);
plot(t,En);
xlabel('$t$','Interpreter','latex');
ylabel('Discrete energy $E^n$','Interpreter','latex');
title(sprintf('Discrete energy, RK4 $\\tau = %.2f$',tau),'Interpreter','latex');

function out = f(c)
out = c - c^3;
end

function Cnext = rk4_step(Cn, tau)
k1 = f(Cn);
k2 = f(Cn + 0.5 * tau * k1);
k3 = f(Cn + 0.5 * tau * k2);
k4 = f(Cn + tau * k3);
Cnext = Cn + tau / 6.0 * (k1 + 2*k2 + 2*k3 + k4);
end

function out = exact_sol(t, c0)
out = (c0 * exp(t)) ./ sqrt(1 - c0^2 + c0^2 * exp(2*t));
end

function out = discrete_energy(C)
out = 1.5 * (C.^2 - 1.0).^2;
end

function [t_arr, C_num, E_num] = simulate_stage3(c0, T, tau)
Nsteps = floor(T / tau);
t_arr = zeros(Nsteps + 1, 1);
C_num = zeros(Nsteps + 1, 1);
E_num = zeros(Nsteps + 1, 1);
C_num(1) = c0;
t_arr(1) = 0.0;
E_num(1) = discrete_energy(C_num(1));
for n = 1:Nsteps
    C_num(n+1) = rk4_step(C_num(n), tau);
    t_arr(n+1) = n * tau;
    E_num(n+1) = discrete_energy(C_num(n+1));
end
end