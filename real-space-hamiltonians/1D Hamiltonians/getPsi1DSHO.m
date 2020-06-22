function [psi] = getPsi1DSHO(x,n)
%[psi] = getPsi1DSHO(x,n)
%Returns analytic solution for 1D SHO wavefunctions. 
run('constants.m');
omega = 2*pi*3e3;
m=mLi;
psi=(2^n*factorial(n))^(-1/2)*(m*omega/(pi*hbar))^(1/4)*exp(-m*omega*x.^2/(2*hbar)).*hermiteH(n,sqrt(m*omega/hbar)*x);
end

