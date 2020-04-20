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
%plot(x,p,':',x,rho,':',x,u,':','LineWidth',2)
%hold on

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
    QL = zeros(3,N);
    QR = zeros(3,N);
    s = zeros(3,N);

    for k = 2:N-1
        for j = 1:3
            s(j,k) = minmod((Qold(j,k+1) - Qold(j,k))/dx,(Qold(j,k) - Qold(j,k-1))/dx);
        end
    end
        %storing QL(i+0.5) at QL(i) and QR(i-0.5) at QR(i-1)
    for k = 1:N-1    
        QL(:,k) = Qold(:,k) + s(:,k).*(0.5*dx);
    end
    for k = 2:N
        QR(:,k-1) = Qold(:,k) - s(:,k).*(0.5*dx);
    end

    %2B) Construct Flux at (i+0.5) and (i-0.5)
    rhoL(:) = QL(1,:);
    rhoR(:) = QR(1,:);
    uL(1,:) = QL(2,:)./rho(1,:); %u1 is saving values of u(i-0.5) in u1(1,:) and u(i+0.5) in u1(2,:) 
    uR(1,:) = QR(2,:)./rho(1,:);
    ETL(:) = QL(3,:);
    ETR(:) = QR(3,:);
    pL(:) = (ETL(:) - ((1/2).*rhoL(:).*(uL(:).^2)))*(gamma - 1);
    pR(:) = (ETR(:) - ((1/2).*rhoR(:).*(uR(:).^2)))*(gamma - 1);

    EL(1,:) = rho(1,:).*uL(1,:);
    ER(1,:) = rho(1,:).*uR(1,:);
    EL(2,:) = EL(1,:).*uL(1,:) + pL(1,:);
    ER(2,:) = ER(1,:).*uR(1,:) + pR(1,:);
    EL(3,:) = ET(1,:).*uL(1,:) + pL(1,:).*uL(1,:);
    ER(3,:) = ET(1,:).*uR(1,:) + pR(1,:).*uR(1,:);

    flux = zeros(3,N); %flux1 is F(i-0.5) and flux2 is F(i+0.5)
    for k = 1:N-1
        c = sqrt(gamma*p(k)./rho(k));
        A = [abs(Qold(2,k)) abs(Qold(2,k)+c) abs(Qold(2,k)-c)];
        al = max(A);
        for j=1:3
            flux(j,k) = 0.5*(ER(j,k) + EL(j,k)) - 0.5*(al)*(QR(j,k) - QL(j,k));
        end
    end


    %2C) Calculate Q at n+1 time
    for k = 2:N-1
        Q(:,k) = Qold(:,k) - C*(flux(:,k) - flux(:,k-1));
    end
    
    %3) Update solution (rho,u,p)
    rho(:) = Q(1,:);
    u(1,:) = Q(2,:)./rho(1,:);
    ET(:) = Q(3,:);
    p(1,:) = (ET(1,:) - ((1/2).*rho(1,:).*(u(1,:).^2)))*(gamma - 1);
    
end

%plot(x,p,x,rho,x,u,'LineWidth',1.5)
plot(x,u,'LineWidth',1.5)
grid on
xlabel('x')
%ylabel('Pressure (Pa), Density (kg/m^3), Velocity (m/s)')
ylabel('Velocity (m/s)')
%legend({'Initial Pressure','Initial Density','Initial Velocity','Final Pressure','Final Density','Final Velocity'},'Fontsize',8,'Location','west')
%plotting analytical solution
%time = 0.174;
%data = analytic_sod(time);
%plot(data.x,data.P,data.x,data.rho,data.x,data.u,'LineWidth',2)
%legend({'Initial Pressure','Initial Density','Initial Velocity','Final Pressure','Final Density','Final Velocity','Analytical Pressure','Analytical Density','Analytical Velocity'},'Fontsize',8,'Location','west')
legend({'Analytical Velocity','Final Velocity'},'Fontsize',12,'Location','west')