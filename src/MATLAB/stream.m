A = 1;
y_max = 10;
x_max = 10;
limits = [1,10];

psi = @(x,y) A*sin(x * pi/ x_max)*sin(y * pi / y_max);

