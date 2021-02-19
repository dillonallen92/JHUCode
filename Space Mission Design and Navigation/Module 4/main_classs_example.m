%initial conditions
muE=3.986e5;       %(km^3/s^2)
%x0Circ=[42241.0800678832;0;0;7.27220521664305e-05];  % (km,rad,km/s,rad/s)
x0Ellipt=[21120.5400339416;0;0;0.000251916578365864]; % (km,rad,km/s,rad/s)


tspan=linspace(0,36*3600,100);

%set up integrator options
integOptions=odeset('abstol',1e-10,'reltol',1e-10);
[t,x]=ode113(@(t,x)two_body_polar(t,x,muE),tspan,x0Ellipt,integOptions);


%plot
figure(1);
hold on;
plot(x(:,1).*cos(x(:,2)),x(:,1).*sin(x(:,2)));
axis equal;
ax=gca;
ax.XLabel.String='x (km)';
ax.YLabel.String='y (km)';



