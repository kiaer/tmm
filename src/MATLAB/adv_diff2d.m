3% Matlab implementation
%
% Christian Ki�r and Anton Almgren

D       = 10; %diffusion constant [m^2 /d]
zmax    = 10;  %depth of water column [m]
zGrid   = 10;  %number of grid cells
dz      = zmax/zGrid; %grid size
z       = linspace(0,zmax,dz);% depth vector located in the middle of each grid cell
xmax    = 10;
xGrid   = 10;
dx      = xmax/xGrid;
x       = linspace(0,xmax,dx);
u_vec   = zeros(zGrid,xGrid);
w_vec   = zeros(zGrid,xGrid);
A       = 1;

syms x y
psi = @(x,y) A*sin(sym(pi) * (x / (xmax - 1))) * sin(sym(pi) * (y / (zmax - 1)));
u = matlabFunction( diff(psi(x,y), y) );
v = matlabFunction( diff(psi(x,y), x) );

for i=1:zGrid
    for j=1:xGrid
        u_vec(i,j) = u(j - 1,i - 1);
        w_vec(i, j) = v(j - 1,i - 1);
    end
end
u_vec(:,end) = 0
w_vec(end,:) = 0

figure
surface(u_vec)
shading interp
figure
surface(w_vec)
shading interp

