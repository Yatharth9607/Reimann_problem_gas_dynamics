%Godunov with MUSCL
clear 
L = 1;
dx = 0.001;
x = (0:dx:L);
N = (L/dx)+1;
T = 0.1644;
dt = 0.0001;
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
plot(x,p,':',x,rho,':',x,u,':')
hold on

%Time steps
for i = 1:n

    %1) Construct Q and E vectors

    Q(1,:) = rho(:);
    Q(2,:) = rho(:).*u(:);
    Q(3,:) = ET(:);

    E(1,:) = rho(:).*u(:);
    E(2,:) = E(1,:).*u(1,:) + p(1,:);
    E(3,:) = ET(:).*u(:) + p(:).*u(:);

    %2) Calculate Q(n+1) using MUSCL scheme
    %2A) Construct u(i+0.5) and u(i-0.5)
    Qold = Q;
    Q1 = zeros(2,N);
    Q2 = zeros(2,N);
    Q3 = zeros(2,N);
    for k = 2:N-1
        for j = 1:3
            s(j,k) = minmod((Qold(j,k+1) - Qold(j,k))/dx,(Qold(j,k) - Qold(j,k-1))/dx);
        end
        Q1(1,k) = Qold(1,k) - s(1,k).*(0.5*dx); 
        Q1(2,k) = Qold(1,k) + s(1,k).*(0.5*dx);
        Q2(1,k) = Qold(2,k) - s(2,k).*(0.5*dx); %Q2 is having values of Q(2,k-0.5) at Q2(1,k) and Q(2,k+0.5) at Q2(2,k)
        Q2(2,k) = Qold(2,k) + s(2,k).*(0.5*dx);
        Q3(1,k) = Qold(3,k) - s(3,k).*(0.5*dx); 
        Q3(2,k) = Qold(3,k) + s(3,k).*(0.5*dx);
    end
    
    %2B) Construct Flux at (i+0.5) and (i-0.5)
%   rho(:) = Q(1,:);
    u1(1,:) = Q2(1,:)./rho(1,:); %u1 is saving values of u(i-0.5) in u1(1,:) and u(i+0.5) in u1(2,:) 
    u1(2,:) = Q2(2,:)./rho(1,:);
%   ET(:) = Q(3,:);
    p1(1,:) = (ET(1,:) - ((1/2).*rho(1,:).*(u1(1,:).^2)))*(gamma - 1);
    p1(2,:) = (ET(1,:) - ((1/2).*rho(1,:).*(u1(2,:).^2)))*(gamma - 1);
    
    Q1(1,:) = rho(1,:);
    Q1(2,:) = rho(1,:);
    Q2(1,:) = rho(1,:).*u1(1,:);
    Q2(2,:) = rho(1,:).*u1(2,:);
    Q3(1,:) = ET(1,:);
    Q3(2,:) = ET(1,:);

    E1(1,:) = rho(1,:).*u1(1,:);
    E1(2,:) = rho(1,:).*u1(2,:);
    E2(1,:) = E1(1,:).*u1(1,:) + p1(1,:);
    E2(2,:) = E1(2,:).*u1(2,:) + p1(2,:);
    E3(1,:) = ET(1,:).*u1(1,:) + p1(1,:).*u1(1,:);
    E3(2,:) = ET(1,:).*u1(2,:) + p1(2,:).*u1(2,:);

    flux1 = zeros(3,N); %flux1 is F(i-0.5) and flux2 is F(i+0.5)
    flux2 = zeros(3,N);
    for k = 1:N-1
        c = sqrt(gamma*p(k)./rho(k));
        A = [abs(Q1(1,k)) abs(Q1(1,k)+c) abs(Q1(1,k)-c)];
        al = max(A);
        flux1(1,k) = 0.5*(E1(1,k) + E1(1,k+1)) - 0.5*(al)*(Q1(1,k+1) - Q1(1,k));
        flux2(1,k) = 0.5*(E1(2,k) + E1(2,k+1)) - 0.5*(al)*(Q1(2,k+1) - Q1(2,k));

        flux1(2,k) = 0.5*(E2(1,k) + E2(1,k+1)) - 0.5*(al)*(Q2(1,k+1) - Q2(1,k));
        flux2(2,k) = 0.5*(E2(2,k) + E2(2,k+1)) - 0.5*(al)*(Q2(2,k+1) - Q2(2,k));

        flux1(3,k) = 0.5*(E3(1,k) + E3(1,k+1)) - 0.5*(al)*(Q3(1,k+1) - Q3(1,k));
        flux2(3,k) = 0.5*(E3(2,k) + E3(2,k+1)) - 0.5*(al)*(Q3(2,k+1) - Q3(2,k));
        
    end

    %2C) Calculate Q at n+1 time
    for k = 2:N-1
        Q(:,k) = Qold(:,k) - C*(flux1(:,k) - flux1(:,k-1));
    end
    
    %{
    ROUGH WORK
    u(x) = u(i) + ((x - x(i))/(x(i+1) - x(i)))*(u(i+1) - u(i));
    du(i)/dt + (1/dx(i))*(F(u(i+0.5)) - F(u(i-0.5))) = 0;
    u(i+0.5) = 0.5*(u(i) + u(i+1));
    u(i-0.5) = 0.5*(u(i-1) + u(i));
    slope limiter:
    s(i) = minmod((u(i+1) - u(i))/dx,(u(i) - u(i-1))/dx);
    minmod(a,b) = sign(a)*min(abs(a),abs(b)) (if a*b > 0)
                  0                          (if a*b <=0)
    u(x) = u(i) + s(i)*(x - x(i));
    %}
    
    %3) Update solution (rho,u,p)
    rho(:) = Q(1,:);
    u(1,:) = Q(2,:)./rho(1,:);
    ET(:) = Q(3,:);
    p(1,:) = (ET(1,:) - ((1/2).*rho(1,:).*(u(1,:).^2)))*(gamma - 1);
    
end

plot(x,p,x,rho,x,u)
legend({'Energy','Density','Velocity'},'Location','southwest')

%{
title(['Fourth Order Runge-Kutta method (C = ',num2str(C),')'])
grid on
xlabel('x')
ylabel('u')
legend({'t = 0','t = 2','t = 4','t = 6'},'Location','southwest')
hold off
%}