%RK with TVD
clear
L = 1;
dx = 0.01;
x = (0:dx:L);
N = (L/dx)+1;
T = 0.1644;
dt = 0.0004;
n = (T/dt);
C = dt/dx;
gamma = 1.4;

%Initial conditions
for i = 1:N
    if x(i) <= 0.5
        rho(i) = 1;
        u(i) = 0;
        p(i) = 2.5*(gamma - 1);
    else
        rho(i) = 0.125;
        u(i) = 0;
        p(i) = 0.25*(gamma - 1);
    end
end
ET(:) = ((1/2).*rho(:).*(u(:).^2)) + (p(:)/(gamma - 1));
%figure(1)
%plot(x,ET,x,rho,x,u)
%legend({'Energy','Density','Velocity'},'Location','southwest')

%Fourth Order RK method
for i = 1:n

    %1) Construct Q and E vectors

    Q(1,:) = rho(:);
    Q(2,:) = rho(:).*u(:);
    Q(3,:) = ET(:);

    E(1,:) = rho(:).*u(:);
    E(2,:) = E(1,:).*u(1,:) + p(1,:);
    E(3,:) = ET(:).*u(:) + p(:).*u(:);

    %2) Calculate Q(n+1) using RK method
    Q1 = Q;
    for j = 1:3
        Ex1(j,:) = Exi(E(j,:),N,dx);
        for l = 1:4
            Ex(j,:) = Exi(E(j,:),N,dx);
            Q(j,:) = Q1(j,:) - (1/(5-l))*dt*(Ex(j,:));
        end
%       Q(j,:) = Q(j,:) + Damp(Q1(j,:),N,0,0.1);
        Q(j,:) = Q(j,:) - C*(sif1(Q1(j,:),N,dx,dt,Ex1(j,:)) - sif2(Q1(j,:),N,dx,dt,Ex1(j,:)));
    end

    %3) Update solution (rho,u,p)
    rho(:) = Q(1,:);
    u(1,:) = Q(2,:)./rho(1,:);
    ET(:) = Q(3,:);
    p(1,:) = (ET(1,:) - ((1/2).*rho(1,:).*(u(1,:).^2)))*(gamma - 1);
    
end

figure(2)
plot(x,ET,x,rho,x,u)
legend({'Energy','Density','Velocity'},'Location','southwest')

%{
title(['Fourth Order Runge-Kutta method (C = ',num2str(C),')'])
grid on
xlabel('x')
ylabel('u')
legend({'t = 0','t = 2','t = 4','t = 6'},'Location','southwest')
hold off
%}