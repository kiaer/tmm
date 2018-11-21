function dydt = ode_npz(t, y, param, i, surf_ind, Ybox)

N = y(1:52749); 
P = y(52749+1:52749*2);
Z = y(52749*2+1:52749*3);

T = i / 2;
%N = matrixToGrid(N, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
%P = matrixToGrid(P, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
%Z = matrixToGrid(Z, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');

g(1:surf_ind) = exp(-(0.025) * param.M).*(1-0.8*sin(pi.*Ybox/180)*cos(2*pi*T/365));



dNdt(1:surf_ind) = - (((sign(N(1:surf_ind)) + 1) ./ 2) .* g'.*( N(1:surf_ind) ./( param.Hp + N(1:surf_ind) )) - ((sign(P(1:surf_ind)) + 1) ./ 2) .* param.r) .* P(1:surf_ind);
dNdt(surf_ind+1:length(N)) =  ((sign(P(surf_ind+1:end)) + 1) ./ 2) .* param.r .* P(surf_ind+1:end);

dPdt(1:surf_ind) = (((sign(N(1:surf_ind)) + 1) ./ 2) .* g'.*( N(1:surf_ind) ./( param.Hp + N(1:surf_ind) )) - ((sign(P(1:surf_ind)) + 1) ./ 2) .* param.r ) .* P(1:surf_ind)...
    - max(0,  ((sign(Z(1:surf_ind)) + 1) ./ 2) .* (param.c .* Z(1:surf_ind) .* (P(1:surf_ind) - param.P0 ) ./ (param.Hz + P(1:surf_ind) - param.P0)));

dPdt(surf_ind+1:length(N)) = - max(0,  ((sign(Z(surf_ind+1:length(N))) + 1) ./ 2) .* (param.c .* Z(surf_ind+1:end) .* (P(surf_ind+1:end) - param.P0 ) ./ (param.Hz + P(surf_ind+1:end) - param.P0)))...
                       - P(surf_ind+1:end) .* ((sign(P(surf_ind+1:end)) + 1) ./ 2) .* param.r;

dZdt(1:surf_ind) = param.eps .* max(0, ((sign(Z(1:surf_ind)) + 1) ./ 2) .* param.c .* Z(1:surf_ind) .*( P(1:surf_ind) - param.P0 ) ./ (param.Hz + P(1:surf_ind) - param.P0)) - ((sign(Z(1:surf_ind)) + 1) ./ 2) .*  param.d .* Z(1:surf_ind);
dZdt(surf_ind+1:length(N)) = param.eps .* max(0, (((sign(Z(surf_ind+1:length(N))) + 1) ./ 2) .* param.c .* Z(surf_ind+1:end) .*( P(surf_ind+1:end) - param.P0 ) ./ (param.Hz + P(surf_ind+1:end)...
    - param.P0))) - ((sign(Z(surf_ind+1:end)) + 1) ./ 2) .*  param.d .* Z(surf_ind+1:end);



%dNdt = gridToMatrix(N, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
%dPdt = gridToMatrix(P, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
%dZdt = gridToMatrix(Z, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
% if T == 100
%     keyboard
% end
dydt = [dNdt, dPdt, dZdt]';

end

