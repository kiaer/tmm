% Calculating incident clear sky solar radiation as a function of time t

% (days) and latitude z (degrees). Relevant constants include:

delto = @(t) -22.44*pi/180*cos(2*pi*(t + 10)/365);

A = @(t,z) cos(delto(t)).*cos(z*pi/180);

B = @(t,z) sin(delto(t)).*sin(z*pi/180);

% to give incident solar radiation So (in units of the average = 340 W / m2) 

So = @(t,z) max(cos(2*pi*t).*A(t,z) + B(t,z),0);

% for which maximum and minimum are

Smax = @(t,z) max(A(t,z) + B(t,z),0);

Smin = @(t,z) max(-A(t,z) + B(t,z),0);


y = linspace(0,365,30000);

figure(1)

clf

plot(y,So(y,90),'-');

