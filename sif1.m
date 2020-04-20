%Roe-Sweby Upward TVD limiter calculator
function [si] = sif1(u,N,dx,dt,Ex)

for i = 2:N-1    
    if 0.5*(u(i+1) + u(i)) == 0
        alph(i) = 0.5*(u(i+1) + u(i));
    else
        alph(i) = (Ex(i+1) - Ex(i))/0.5*(u(i+1) + u(i));
    end
    a = sign(alph(i));
    r(i) = (u(i+1+a) - u(i+a))/(0.5*(u(i+1) + u(i)));
    G(i) = (sign(1) + sign(r(i))/2)*(min(abs(1),abs(r(i))));
    si(i) = (0.5*G(i)*(abs(alph(i)) + (dt/dx)*(alph(i).^2)) - abs(alph(i)))*(0.5*(u(i+1) + u(i)));
end
si(N) = 0;