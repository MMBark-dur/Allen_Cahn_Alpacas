epsilon = 0.01;
a = -2;
b = 4;
T = 10;

h = 0.01; % space step (adjustable)
dt = 0.01; % time step (adjustable)
%dt = 0.005 
%dt=0.0025
x = linspace(a, b, round((b-a)/h) + 1);

u0 = 0.5 * (1 - tanh(x / (2 * sqrt(2) * epsilon)));
u = u0;

for s = 1:round(T/dt)   %s: time step index
    u1 = u;

    for g = 2:length(x)-1   %g: spatial grid-point index
        uxx = (u(g-1) - 2*u(g) + u(g+1)) / h^2;
        rate = epsilon^2 * uxx + u(g) - u(g)^3;
        u1(g) = u(g) + dt * rate;
    end

    u = u1;
end

%check plot
figure(1)
plot(x, u0, '--', x, u, '-')
xlabel('x')
ylabel('u')
legend('Initial profile', 'Final profile')
grid on


%numerical wave speed c
j = find(u(1:end-1) > 0.5 & u(2:end) <= 0.5, 1);
xmid = x(j) + (0.5 - u(j)) * (x(j+1) - x(j)) / (u(j+1) - u(j));
%xmid: final wave centre
cnum = xmid / T;

fprintf('Final wave centre: %.6f\n', xmid);
fprintf('Numerical wave speed: %.6f\n', cnum);

%Final wave centre: 0.212126  Numerical wave speed: 0.021213 (This is it
%for now)

cana = 3 * sqrt(2) * epsilon / 2; % analytic wave speed
uref = 0.5 * (1 - tanh((x - cana * T) / (2 * sqrt(2) * epsilon)));
%reference profile at t=10

figure(2)
plot(x, u, '-', x, uref, '--')
xlabel('x')
ylabel('u')
legend('Forward Euler', 'Analytic reference')
title('Comparison at t = 10')
grid on

% maximum absolute error
errmax = max(abs(u - uref));
fprintf('Maximum absolute error at t = %.1f: %.6e\n', T, errmax);