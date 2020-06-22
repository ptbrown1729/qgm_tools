function [H] = rabiH(t)
%RABIH Summary of this function goes here
%   Detailed explanation goes here
run('constants.m');

omega=omega_D2;
rabiFrq=2*pi*3e3;
Eexcited=omega_D2;%+2*pi*1e6;
H=hbar*[[0,rabiFrq*cos(omega*t)];[rabiFrq*cos(omega*t),Eexcited]];

%expect splitting to be sqrt(Delta^2+4*rabiFrq^2);

% state_i=[1 0];
% state_t=@(t) state_i*V(:,1)*V(:,1)*exp(1i*D(1,1)*t)+state_i*V(:,2)*V(:,2)*exp(1i*D(2,2)*t);
% trans_prob=@(t) [0 1]*state_t(t);
end

