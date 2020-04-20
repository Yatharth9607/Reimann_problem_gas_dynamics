clear
a = 0;
for N = 100:100:1000
    a = a + 1;
    b(a) = N;
    L = 1;
    C = 0.1;
    T = 0.1644;
    gamma = 1.4;
    dx = L/N;
    x = linspace(0,L,N);
    dt = C*dx;
    n = (T/dt);
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
        ET(i) = ((1/2).*rho(i).*(u(i).^2)) + (p(i)/(gamma - 1));
    end

    %Fourth Order RK method
    tic
    for i = 1:n

        %1) Construct Q and E vectors

        Q = zeros(3,N);        
        Q(1,:) = rho(:);
        Q(2,:) = rho(:).*u(:);
        Q(3,:) = ET(:);
        
        E = zeros(3,N);
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

    t(a) = toc;
%{
    plot(x,u,'LineWidth',1)
    grid on
    xlabel('x (m)')
    ylabel('Velocity (m/s)')
    hold on
%}
end
%hold off
%legend({'dx = 0.01 m','dx = 0.002 m','dx = 0.001 m'},'Fontsize',8,'Location','west')
figure(1)
plot(b,t)
hold on
%plotting analytical solution
%time = 0.174;
%data = analytic_sod(time);
%plot(data.x,data.P,data.x,data.rho,data.x,data.u,'LineWidth',2)
%legend({'Initial Pressure','Initial Density','Initial Velocity','Final Pressure','Final Density','Final Velocity','Analytical Pressure','Analytical Density','Analytical Velocity'},'Fontsize',8)