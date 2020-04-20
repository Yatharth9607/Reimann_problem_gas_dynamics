clear
L = 1;
dx = 0.01;
x = (0:dx:L);
N = (L/dx)+1;
T = 0.16;
dt = 0.001;
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
figure(1)
plot(x,ET,x,rho,x,u)
legend({'Energy','Density','Velocity'},'Location','southwest')
%{%
%McCormac scheme
for i = 1:20
    %predictor step
    %1) Construct Q and E vectors

    Q(1,:) = rho(:);
    Q(2,:) = rho(:).*u(:);
    Q(3,:) = ET(:);

    E(1,:) = rho(:).*u(:);
    E(2,:) = E(1,:).*u(1,:) + p(1,:);
    E(3,:) = ET(:).*u(:) + p(:).*u(:);

    %2) Calculate Q* using upwinding scheme
    Q1 = Q;
    for j = 1:N-1
        Q(:,j) = Q1(:,j) - C*(E(:,j+1) - E(:,j));
    end

    %3) Update solution (rho,u,p,ET)
    rho(:) = Q(1,:);
    u(1,:) = Q(2,:)./Q(1,:);
    ET(:) = Q(3,:);
    p(1,:) = (ET(1,:) - ((1/2).*rho(1,:).*(u(1,:).^2)))*(gamma - 1);

    %corrector step
    %1) Construct E vector
    E(1,:) = rho(:).*u(:);
    E(2,:) = E(1,:).*u(1,:) + p(1,:);
    E(3,:) = ET(:).*u(:) + p(:).*u(:);

    %2) Calculate Q(n+1) using upwinding scheme
    for j = 2:N
        Q(:,j) = (1/2)*(Q1(:,j) + Q(:,j) - C*(E(:,j) - E(:,j-1)));
    end
    
    %3) Update solution (rho,u,p)
    rho(:) = Q(1,:);
    u(:) = Q(2,:)./Q(1,:);
    ET(:) = Q(3,:);
    p(1,:) = (ET(1,:) - ((1/2).*rho(1,:).*(u(1,:).^2)))*(gamma - 1);
    
end
%}

figure(2)
plot(x,ET,x,rho,x,u)
legend({'Energy','Density','Velocity'},'Location','southwest')