%Ex calculator
function [Ex] = Exi(u1,N,dx)

for i = 2:N-1
    Ex(i) = ((u1(i+1).^2) - (u1(i-1).^2))/(4*dx);
end
Ex(N) = 0;