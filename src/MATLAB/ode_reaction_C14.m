function [dydt] = ode_reaction_C14(t,Y,param)
%ODE_ADV_DIFF2D Summary of this function goes here
%   Detailed explanation goes here

%H = 10; %half life 5730 yrs

D = 1/(2*param.H); %decay const (decay pr year)

C = reshape(Y, param.zGrid, param.xGrid);


%advective flux 
for i = 1:param.zGrid
    for j = 2:param.xGrid
        if param.u_vec(i,j) >= 0
            Jax(i,j) = (param.u_vec(i,j-1) * C(i,j-1))*param.dz;
           Jax(i,1) = 0;
        else
            Jax(i,j) = (param.u_vec(i,j) * C(i,j))*param.dz;
            Jax(i,param.xGrid+1) = 0;
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
            Jaz(param.zGrid+1,j) = 0;
        end
        %Jdz(2:param.zGrid,j) = -param.D*((C(2:param.zGrid,j)-C(1:param.zGrid-1,j))./param.dz)*param.dx ;
        
    end
end

for j = 1:param.xGrid
    Jdz(2:param.zGrid,j) = -param.D*((C(2:param.zGrid,j)-C(1:param.zGrid-1,j))./param.dz)*param.dx ;
end

% Jaz([1 param.zGrid+1],:) = 0;
% Jax(:,[1 param.xGrid+1]) = 0;
 Jdz(param.zGrid + 1,:) = 0;
 Jdz(1,:) =  -param.D .* ((Jdz(1) - param.atm) ./ param.dz);
% Jaz(1,:) = 0;
% Jax(:,1) = 0;
% Jdz(1,:) = 0;

Jx = Jax; %+ Jdx
Jz = Jaz + Jdz;
dcdt = (Jx(:,1:param.xGrid)-Jx(:,2:param.xGrid+1))/param.dx ...
    + (Jz(1:param.zGrid,:)-Jz(2:param.zGrid+1,:))/param.dz...
    - D * C;

dydt = dcdt(:);
end

