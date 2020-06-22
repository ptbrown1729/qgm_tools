%test dampedOsc1D, supposed to solve dt^2n + (2/l) * dtn + omega^2 * n = 0
l = 3;
omega_o = 6;
xo = 2;
vo = 0;

%function
y_fn = @(t) dampedOsc1D([xo, vo, omega_o, l, 0, 0], t);


%numerical solution
dy = @(t,y) [y(2); -y(2) * 2/l - omega_o^2 * y(1)];
tmax = 2*pi/omega_o*10;
[ts,ys]=ode45(dy,[0,tmax],[xo, vo]);

figure;
plot(ts, ys(:,1),'b.');
hold on;
interp_t = linspace(min(ts), max(ts), 1000);
plot(interp_t, y_fn(interp_t),'r');
plot(interp_t, exp1D([0, xo, l, 0], interp_t));
% plot(interp_t, exp1D([0, xo, omega_o^2*l/2, 0], interp_t)); 
grid on;
xlabel('t');
ylabel('y');
legend({'numerical solution','analytic soln', 'decaying expt @ l'});
ttl = sprintf('Solution to dt^2*n + (2/l) dt*n + omega_o^2*n = 0 \n l = %0.2f, wo = %0.2f \n xo = %0.2f, vo = %0.2f', l, omega_o, xo, vo);
suptitle(ttl);