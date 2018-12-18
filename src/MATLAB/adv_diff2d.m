%% Matlab implementation
%
% Christian Kiaer and Anton Almgren
%%
clc
clear all
close all
%%

param.D       = 0.0005; %diffusion constant [km^2 /yr]
param.zmax    = 4;  %depth of water column [m]
param.zGrid   = 20;  %number of grid cells

param.dz      = param.zmax/param.zGrid;
param.z       = [0.5*param.dz:param.dz:param.zmax]';% depth vector located in the middle of each grid cell

param.xmax    = 5000;
param.xGrid   = 20;

param.dx      = param.xmax/param.xGrid;
param.x       = [0.5*param.dx:param.dx:param.xmax]';

param.u_vec   = zeros(param.zGrid,param.xGrid);
param.w_vec   = zeros(param.zGrid,param.xGrid);

param.atm     = 1;   %Atm concentration
param.H       = 10;  %Half

A             = 1000;



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
xlabel('South - North [km]')
ylabel('depth [km]')
c=colorbar('FontSize',14);
c.Label.String='velocity [km/yr]';
set(gca,'FontSize', 14)
print('../../fig/2d_U', '-dpng', '-r300')


figure
surface(param.x,param.z,param.w_vec)
axis ij
shading interp
title('w')
c=colorbar('FontSize',14);
xlabel('South - North [km]')
ylabel('depth [km]')
c.Label.String='velocity [km/yr]';
set(gca,'FontSize', 14)
print('../../fig/2d_W', '-dpng', '-r300')


%%
tspan = 0:1;
C0 = zeros(param.zGrid,param.xGrid);
%C0(2*param.zGrid/10,5*param.xGrid/10) = 10;
C0(1,:) = 1;


%C0(2,end-1) = 10;
tflag = tic;
[t,Y] = ode23tb(@ode_adv_diff2d, [0 10], C0(:), [], param);
Solvetime_ODE=toc(tflag)

C = reshape(Y',param.zGrid,param.xGrid,length(t));



figure
surface(param.x,param.z,C(:,:,end))
axis ij
shading flat
colorbar
title(['ODE solution 10 yrs, Sol. time =',num2str(Solvetime_ODE),'s'])
ylabel('Depth [km]')
xlabel('South - North [km]')
caxis([0 0.8])
set(gca,'FontSize', 16)

%% Transport Matrix
%load('../../bin/A_20x20.mat')
A = sparse(param.xGrid*param.zGrid, param.xGrid*param.zGrid);
ind = 1;
tflag = tic;
for i=1:param.zGrid
    for j=1:param.xGrid
        Ca0 = zeros(param.zGrid, param.xGrid);
        Ca0(j,i) = 1;
        %Ca0
        [t,Y] = ode23tb(@ode_adv_diff2d, [0 0.1], Ca0(:), [], param);
        
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
Buildingtime_A =toc(tflag); %~200s
%%

C0 = zeros(param.zGrid,param.xGrid);
C0(1,:) = 1;
An = A^10;

tflag = tic;
Ca = An^10 * C0(:);
Ca = reshape(Ca, param.zGrid, param.xGrid);
Solvetime_TMM = toc(tflag);


figure
surface(param.x,param.z, Ca)
axis ij
shading flat
sum(sum(Ca))
colorbar
xlabel('South - North [km]')
ylabel('depth [km]')
title(['TMM solution 10 yrs, Sol. time =',num2str(Solvetime_TMM),'s'])
caxis([0 0.8])
set(gca,'FontSize', 16)


%%
An = A ^ 10;
C0 = zeros(param.zGrid,param.xGrid);
%C0(2*param.zGrid/10,5*param.xGrid/10) = 10;
C0(1,:) = 1;
Cn = C0;
D = 1/(2*param.H);
tflag = tic;
for i = 1:100
    %Cn = Cn(:)
    Cn(1,:) = param.atm;
    Cn = An * Cn(:) - D .* Cn(:);
    %size(temp)
    Cn = reshape(Cn, param.zGrid, param.xGrid);
    %Cn(1,:) = -param.D .* ((Cn(1,:) - param.atm) ./ param.dz) .* param.dx;
end
Solvetime_TMMreaction = toc(tflag);


figure
contourf(param.x,param.z,Cn,5)
axis ij
c=colorbar;
c.Label.String='concentration';
caxis([0 1.2])
shading flat
xlabel('South - North [km]')
ylabel('depth [km]')
title(['TMM solution 10 yrs, Sol. time =',num2str(Solvetime_TMMreaction),'s'])
set(gca,'FontSize', 14)
print('../../fig/2d_tmm_sol', '-dpng', '-r300')

%%
tspan = 0:100;
C0 = zeros(param.zGrid,param.xGrid);
C0(1,:) = 1;

tflag = tic;
[t,Y] = ode23tb(@ode_reaction_C14, tspan, C0(:), [], param);
C = reshape(Y',param.zGrid,param.xGrid,length(t));
Solvetime_ODEreaction = toc(tflag);

figure
contourf(param.x,param.z,C(:,:,end),5)
axis ij
shading flat
c=colorbar;
caxis([0 1.2])
xlabel('South - North [km]')
ylabel('depth [km]')
c.Label.String='concentration';
title(['ODE solution 10 yrs, Sol. time =',num2str(Solvetime_ODEreaction),'s'])
set(gca,'FontSize', 14)
print('../../fig/2d_ode_sol', '-dpng', '-r300')


figure
contourf(param.x,param.z,C(:,:,end)-Cn,5)
axis ij
shading flat
c=colorbar;
%caxis([0 1.2])
xlabel('South - North [km]')
ylabel('depth [km]')
c.Label.String='difference in concentration';
title('ODE solution - TMMsolution')
set(gca,'FontSize', 14)
print('../../fig/2d_ode_tm_diff', '-dpng', '-r300')

%%
figure
spy(A)
set(gca,'FontSize', 14)
print('../../fig/sparsity_2d_ideal', '-dpng', '-r300')

%title('Structure of transport matrix')

%%
[X,Z]=meshgrid(param.x,param.z);
figure
streamslice(X,Z,param.u_vec,param.w_vec)
axis ij
%ylim([0 4])
%%
figure
streamline(X,Z,param.u_vec,param.w_vec)
axis ij

