%minmod function
function [M] = minmod(x1,x2)

A = [abs(x1) abs(x2)];

M = 0.5*(sign(x1) + sign(x2))*min(A);