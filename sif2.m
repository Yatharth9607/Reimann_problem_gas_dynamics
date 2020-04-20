%Roe-Sweby Upward TVD limiter calculator
function [si] = sif2(u,N,dx,dt,Ex)

for i = 2:N-1    
    if 0.5*(u(i) + u(i-1)) == 0
        alph(i) = 0.5*(u(i) + u(i-1));
    else
        alph(i) = (Ex(i) - Ex(i-1))/0.5*(u(i) + u(i-1));
    end
    a = sign(alph(i));
    r(i) = (u(i+a) - u(i-1+a))/(0.5*(u(i) + u(i-1)));
    G(i) = (sign(1) + sign(r(i))/2)*(min(abs(1),abs(r(i))));
    si(i) = (0.5*G(i)*(abs(alph(i)) + (dt/dx)*(alph(i).^2)) - abs(alph(i)))*(0.5*(u(i) + u(i-1)));
end
si(N) = 0;