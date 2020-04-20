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

plot(x,p,':',x,rho,':',x,u,':','LineWidth',2)
hold on
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
    Ex = zeros(1,N);
    Q1 = Q;
    for j = 1:3
        for l = 1:4
            for k = 2:N-1
                Ex(j,k) = 0.5*(E(j,k+1) - E(j,k-1))/dx;
            end
            Q(j,:) = Q1(j,:) - (1/(5-l))*dt*(Ex(j,:));
            
            rho(:) = Q(1,:);
            u(1,:) = Q(2,:)./rho(1,:);
            ET(:) = Q(3,:);
            p(1,:) = (ET(1,:) - ((1/2).*rho(1,:).*(u(1,:).^2)))*(gamma - 1);
            
            E(1,:) = rho(:).*u(:);
            E(2,:) = E(1,:).*u(1,:) + p(1,:);
            E(3,:) = ET(:).*u(:) + p(:).*u(:);
        end
        Q(j,:) = Q(j,:) + Damp(Q1(j,:),N,0.05,0);
    end

    %3) Update solution (rho,u,p)
    rho(:) = Q(1,:);
    u(1,:) = Q(2,:)./rho(1,:);
    ET(:) = Q(3,:);
    p(1,:) = (ET(1,:) - ((1/2).*rho(1,:).*(u(1,:).^2)))*(gamma - 1);
    
end

plot(x,p,x,rho,x,u,'LineWidth',1.5)
grid on
xlabel('x')
ylabel('Pressure (Pa), Density (kg/m^3), Velocity (m/s)')
%legend({'Initial Pressure','Initial Density','Initial Velocity','Final Pressure','Final Density','Final Velocity'},'Fontsize',8,'Location','west')
%plotting analytical solution
%time = 0.174;
%data = analytic_sod(time);
%plot(data.x,data.P,data.x,data.rho,data.x,data.u,'LineWidth',2)
legend({'Analytical Pressure','Analytical Density','Analytical Velocity','Initial Pressure','Initial Density','Initial Velocity','Final Pressure','Final Density','Final Velocity'},'Fontsize',8,'Location','west')