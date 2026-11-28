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
    uSIE = u0;
    uFE = u0; % ADDED

    r = epsilon^2*dt/h^2;
    e = ones(N-1,1);
    A = spdiags([-r*e, (1+2*r)*e, -r*e], [-1,0,1], N-1,N-1);

    UFE = zeros(M+1,N+1);
    UFE(1,:) = u0;
    maxErrorFE = zeros(M+1,1);
    maxErrorFE(1) = max(abs(u0-q(x,t(1))));

    USIE = zeros(M+1,N+1);
    USIE(1,:) = u0;
    maxErrorSIE = zeros(M+1,1);
    maxErrorSIE(1) = max(abs(u0-q(x,t(1))));

    for s = 1:M
        u1 = uSIE;
        rhs = (uSIE(2:end-1) + dt*(uSIE(2:end-1)-uSIE(2:end-1).^3))';

        u1(1) = q(a,t(s+1));
        u1(end) = q(b,t(s+1));
        rhs(1) = rhs(1) + r*u1(1);
        rhs(end) = rhs(end) + r*u1(end);

        u1(2:end-1) = (A \ rhs)';
        uSIE = u1;
        USIE(s+1,:) = uSIE;
        maxErrorSIE(s+1) = max(abs(uSIE-q(x,t(s+1))));
    end

    for s = 1:M
        u1 = uFE;
        u1(1) = q(a,t(s+1));
        u1(end) = q(b,t(s+1));

        for g = 2:length(x)-1
            uxx = (uFE(g-1) - 2*uFE(g) + uFE(g+1)) / h^2;
            rate = epsilon^2 * uxx + uFE(g) - uFE(g)^3;
            u1(g) = uFE(g) + dt * rate;
        end

        uFE = u1;
        UFE(s+1,:) = uFE; % ADDED
        maxErrorFE(s+1) = max(abs(uFE-q(x,t(s+1)))); % CORRECTED
    end

    uExact = q(x,T);

    figure(1);
    subplot(2,3,p)
    plot(x,uFE,'-',x,uSIE,'--',x,uExact,':', 'LineWidth',3);
    xlabel('x');
    ylabel('u');
    ylim([-1,1])
    title(sprintf('h=%g, dt=%g',h,dt));
    set(gca, 'fontsize', 20)
    grid on;

    figure(2);
    subplot(2,3,p)
    plot(t,maxErrorFE,'-',t,maxErrorSIE,'--', 'LineWidth',3);
    xlabel('Time');
    ylabel('Max abs err');
    ylim([0,1])
    xlim([0,10])
    title(sprintf('h=%g, dt=%g',h,dt));
    set(gca, 'fontsize', 20)
    grid on;
end

figure(1);
sgtitle(sprintf('Forward Euler, Semi-implicit Euler and exact solution at T = %g',T));

figure(2);
sgtitle('Forward vs Semi-implicit Euler: maximum absolute error over time');