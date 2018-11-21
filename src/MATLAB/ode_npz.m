function dydt = ode_npz(t, y, param, i, surf_ind, Ybox, j, layerind,layer)

% N = y(1); 
% P = y(2);
% Z = y(3);
N = y(1:layer(j));
P = y(layer(j)+1:2*layer(j));
Z = y(2*layer(j)+1:3*layer(j));

T = i / 2;
% N = matrixToGrid(N, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
% P = matrixToGrid(P, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
% Z = matrixToGrid(Z, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');

g(1:layer(1)) = exp(-(0.025) * param.M).*(1-0.8*sin(pi.*Ybox(1:surf_ind)/180)*cos(2*pi*T/365));

% if j ==1
%     g(1:layer(1)) = exp(-(0.025) * param.M).*(1-0.8*sin(pi.*Ybox(1:layer(1))/180)*cos(2*pi*T/365));
% elseif j==2
%     g(1:layer(2)) = exp(-(0.025) * param.M).*(1-0.8*sin(pi.*Ybox(layer(1)+1:surf_ind)/180)*cos(2*pi*T/365));
% end
if j ==1%||j==2
    dNdt = - (((sign(N) + 1) ./ 2) .* g(j) .* ( N ./( param.Hp + N)) - ((sign(P) + 1) ./ 2) .* param.r) .* P;

    dPdt = (((sign(N) + 1) ./ 2) .* g(j).* ( N ./( param.Hp + N )) - ((sign(P) + 1) ./ 2) .* param.r ) .* P ...
    - max(0,  ((sign(Z) + 1) ./ 2) .* (param.c .* Z .* (P - param.P0 ) ./ (param.Hz + P - param.P0)));

    dZdt = param.eps .* max(0, ((sign(Z) + 1) ./ 2) .* param.c .* Z .*( P - param.P0 ) ./ (param.Hz + P - param.P0)) - ((sign(Z) + 1) ./ 2) .*  param.d .* Z;
else
    
    dNdt =  ((sign(P) + 1) ./ 2) .* param.r .* P;


    dPdt = - max(0,  ((sign(Z) + 1) ./ 2) .* (param.c .* Z .* (P - param.P0 ) ./ (param.Hz + P - param.P0)))...
                       - P .* ((sign(P) + 1) ./ 2) .* param.r;

    dZdt = param.eps .* max(0, (((sign(Z) + 1) ./ 2) .* param.c .* Z .*( P - param.P0 ) ./ (param.Hz + P...
                  - param.P0))) - ((sign(Z) + 1) ./ 2) .*  param.d .* Z;

end




% dNdt = gridToMatrix(N, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
% dPdt = gridToMatrix(P, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
% dZdt = gridToMatrix(Z, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
% if T == 100
%     keyboard
% end
dydt = [dNdt; dPdt; dZdt];

end

