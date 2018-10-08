function [dydt] = ode_adv_diff2d(t,Y,param)
%ODE_ADV_DIFF2D Summary of this function goes here
%   Detailed explanation goes here

C = reshape(Y, param.zGrid, param.xGrid);



%advective flux 
for i = 1:param.zGrid
    for j = 2:param.xGrid
        if param.u_vec(i,j) >= 0
            Jax(i,j) = (param.u_vec(i,j-1) * C(i,j-1))*param.dz;
           Jax(i,1) = 0;
        else
            Jax(i,j) = (param.u_vec(i,j) * C(i,j))*param.dz;
            Jax(i,param.zGrid+1) = 0;
            %diff flux(z-dir)
            
        end
        

    end
end
for i = 2:param.zGrid
    for j = 1:param.xGrid
        if param.w_vec(i,j) >= 0 
            Jaz(i,j) = (param.w_vec(i-1,j) * C(i-1,j))*param.dx;
            Jaz(1,j) = 0;
        else
            Jaz(i,j) = (param.w_vec(i,j) * C(i,j))*param.dx;
            Jaz(param.xGrid+1,j) = 0;
        end
        Jdz(2:param.zGrid,j) = -param.D*((C(2:param.zGrid,j)-C(1:param.zGrid-1,j))./param.dz)*param.dx ;
        
    end
end

% Jaz([1 param.zGrid+1],:) = 0;
% Jax(:,[1 param.xGrid+1]) = 0;
 Jdz([1 param.zGrid+1],:) = 0;
% Jaz(1,:) = 0;
% Jax(:,1) = 0;
% Jdz(1,:) = 0;

Jx = Jax; %+ Jdx
%Ix = find(Jax < 0.000000001);
%Jx(Ix) = 0;
Jz = Jaz + Jdz;
%Iz = find(Jax < 0.000000001);
%Jz(Iz) = 0;

dcdt = (Jx(:,1:param.xGrid)-Jx(:,2:param.xGrid+1))/param.dx ...
    + (Jz(1:param.zGrid,:)-Jz(2:param.zGrid+1,:))/param.dz;

dydt = dcdt(:);
end

