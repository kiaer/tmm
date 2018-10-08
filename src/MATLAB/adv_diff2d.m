%% Matlab implementation
%
% Christian Kiaer and Anton Almgren
%%
clc
clear all
close all
%%

param.D       = 1; %diffusion constant [m^2 /d]
param.zmax    = 4;  %depth of water column [m]
param.zGrid   = 40;  %number of grid cells

param.dz      = param.zmax/param.zGrid;
param.z       = [0.5*param.dz:param.dz:param.zmax]';% depth vector located in the middle of each grid cell

param.xmax    = 5000;
param.xGrid   = 60;

param.dx      = param.xmax/param.xGrid;
param.x       = [0.5*param.dx:param.dx:param.xmax]';

param.u_vec   = zeros(param.zGrid,param.xGrid);
param.w_vec   = zeros(param.zGrid,param.xGrid);


A       = 1000;



syms x y
psi = @(x,y) A*sin(sym(pi) * (x^2 / (param.xmax^2))) * sin(sym(pi) * (sqrt(y) / sqrt(param.zmax)));
u = matlabFunction( diff(psi(x,y), y) );
v = matlabFunction( -diff(psi(x,y), x) );

%henter u og w hastighederne i de tilsvarende x- og z- koordinater
for i=1:length(param.z)
    for j=1:length(param.x)
        param.u_vec(i,j) = u(param.x(j),param.z(i));
        param.w_vec(i, j) = v(param.x(j),param.z(i));
    end
end

param.u_vec(1,1) = param.u_vec(1,2);
param.u_vec(end,end) = param.u_vec(end,end-1);
param.w_vec(1,end) = param.w_vec(2,end);
param.w_vec(end,1) = param.w_vec(end-1,1);




% for i=1:param.zGrid
%     for j=1:param.xGrid
%         param.u_vec(i,j) = u(j - 1,i - 1);
%         param.w_vec(i, j) = v(j - 1,i - 1);
%     end
% end
%param.u_vec(:,end) = 0;
%param.w_vec(end,:) = 0;
%%
figure
surface(param.x,param.z,param.u_vec)
shading interp
axis ij
title('u')
colorbar

figure
surface(param.x,param.z,param.w_vec)
axis ij
shading interp
title('w')
colorbar
%%
tspan = 0:1;
C0 = zeros(param.zGrid,param.xGrid);
%C0(2*param.zGrid/10,5*param.xGrid/10) = 10;
C0(1,1) = 1;

%C0(2,end-1) = 10;
tic
[t,Y] = ode23tb(@ode_adv_diff2d, [0 0.001], C0(:), [], param);
toc

C = reshape(Y',param.zGrid,param.xGrid,length(t));
%% Transport Matrix

A = sparse(param.xGrid*param.zGrid, param.xGrid*param.zGrid);
ind = 1;
tic
for i=1:param.zGrid
    for j=1:param.xGrid
        Ca0 = zeros(param.zGrid, param.xGrid);
        Ca0(j,i) = 1;
        %Ca0
        [t,Y] = ode23tb(@ode_adv_diff2d, [0 0.001], Ca0(:), [], param);
        
        if mod(ind,10) == 0
            ind
        end
       % Ynon = nonzeros(Y(end,:))
        Y_e = Y(end,:)';
        %sum_y = sum(Y_e);
        Y_e(Y_e <= 0.005) = 0;
        Y_e = Y_e / sum(Y_e);
        %Y_e = Y_e / sum(Y_e);
        %A = sparse(1:param.xGrid*param.zGrid,ind,Y_e,param.xGrid*param.zGrid, param.xGrid*param.zGrid);
        A(:,ind) = Y_e;
        ind = ind + 1;
        %Ca0(j,i) = 0;
    end
end
toc
%%
tspan = 0:1;
C0 = zeros(param.zGrid,param.xGrid);
%C0(2*param.zGrid/10,5*param.xGrid/10) = 10;
C0(1,:) = 1;

An = A^100 * C0(:)
An = reshape(An, param.zGrid, param.xGrid)
surface(param.x,param.z, An)
axis ij
shading interp
sum(An, 'all')
%%
for i = length(t):-10:1
    figure
    surface(param.x,param.z,C(:,:,i))
    shading interp
    colorbar
    caxis([0 1])
    axis ij
    xlabel('x')
    ylabel('depth')
    title(['t =' num2str(i)])
end
%%
figure
surface(param.x,param.z,C(:,:,end))
shading interp
axis ij
xlabel('x')
ylabel('depth')
    
