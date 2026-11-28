epsilon = 0.01;
a = -2;
b = 4;
T = 10;

hs = [0.005, 0.01, 0.5];
dts = [0.005, 0.01, 0.5];
pairs = [1,3; 2,2; 2,3; 3,1; 3,2; 3,3];

for p = 1:size(pairs,1)
    i = pairs(p,1);
    k = pairs(p,2);

    h = hs(i);
    dt = dts(k);

    N = round((b-a)/h);
    M = round(T/dt);
    x = linspace(a,b,N+1);
    h = (b-a)/N;
    dt = T/M;
    t = linspace(0,T,M+1);
    
    c = 3*sqrt(2)*epsilon/2;
    q = @(x,t) 0.5*(1-tanh((x-c*t)/(2*sqrt(2)*epsilon)));
    u0 = q(x,0);
    u = u0;
    
    r = epsilon^2*dt/h^2;
    e = ones(N-1,1);
    A = spdiags([-r*e, (1+2*r)*e, -r*e], [-1,0,1], N-1,N-1);
    
    U = zeros(M+1,N+1);
    U(1,:) = u0;
    maxError = zeros(M+1,1);
    
    for s = 1:M
        u1 = u;
        rhs = (u(2:end-1) + dt*(u(2:end-1)-u(2:end-1).^3))';
    
        u1(1) = q(a,t(s+1));
        u1(end) = q(b,t(s+1));
        rhs(1) = rhs(1) + r*u1(1);
        rhs(end) = rhs(end) + r*u1(end);
    
        u1(2:end-1) = (A \ rhs)';
        u = u1;
        U(s+1,:) = u;
        maxError(s+1) = max(abs(u-q(x,t(s+1))));
    end
    
    uExact = q(x,T);
    
    figure(1);
    subplot(2,3,p)
    plot(x,u0,':',x,u,'-',x,uExact,'--', 'LineWidth', 3);
    xlabel('x');
    ylabel('u');
    title(sprintf('h=%g, dt=%g',h,dt));
    set(gca, 'fontsize', 20)
    grid on;
    
    figure(2);
    subplot(2,3,p)
    plot(t,maxError,'-', 'LineWidth', 3);
    xlabel('Time');
    ylabel('Max abs err');
    title(sprintf('h=%g, dt=%g',h,dt));
    set(gca, 'fontsize', 20)
    ylim([0,1])
    grid on;
    
    j = find(u(1:end-1) > 0.5 & u(2:end) <= 0.5,1);
    if isempty(j)
        xmid = NaN;
    else
        xmid = x(j) + (0.5-u(j))*(x(j+1)-x(j))/(u(j+1)-u(j));
    end
    cnum = xmid/T;
    
    fprintf('Final wave centre: %.6f\n',xmid);
    fprintf('Numerical wave speed: %.6f\n',cnum);
    fprintf('Maximum absolute error at T: %.8g\n',maxError(end));
end

figure(1);
sgtitle(sprintf('Semi-implicit Euler and exact solution at T = %g',T));
legend('Initial profile','Numerical solution at T','Exact solution at T');

figure(2);
sgtitle('Semi-implicit Euler: maximum absolute error over time');


