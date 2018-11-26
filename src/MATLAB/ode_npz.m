function dydt = ode_npz(t, y, param, i, surf_ind, Ybox, j, layer)

% Split vector into N, P and Z
N = y(1:layer(j));
P = y(layer(j)+1:2*layer(j));
Z = y(2*layer(j)+1:3*layer(j));

% Time in days
T = i / 2;

% If surface layer
if j == 1
    g(1:layer(j)) = exp(-(0.025) * param.M).*(1-0.8*sin(pi.*Ybox(1:surf_ind)/180)*cos(2*pi*T/365));

    dNdt = - (((sign(N) + 1) ./ 2) .* (g' .* ( N ./( param.Hp + N))) - ((sign(P) + 1) ./ 2) .* param.r) .* P;

    dPdt = (((sign(N) + 1) ./ 2) .* g' .* ( N ./( param.Hp + N )) - ((sign(P) + 1) ./ 2) .* param.r ) .* P ...
           - ((sign(Z) + 1) ./ 2) .* (param.c .* Z .* (((sign(P) + 1) ./ 2) .* P - param.P0 ) ./ (param.Hz + ((sign(P) + 1) ./ 2) .* P - param.P0));

    dZdt = param.eps .* ((sign(Z) + 1) ./ 2) .* param.c .* Z .*( ((sign(P) + 1) ./ 2) .* P - param.P0 ) ./ (param.Hz + ((sign(P) + 1) ./ 2) .* P - param.P0) - ((sign(Z) + 1) ./ 2) .*  param.d .* Z;
    
% Every other layer
else
    
    dNdt =  ((sign(P) + 1) ./ 2) .* param.r .* P;


    dPdt = - ((sign(Z) + 1) ./ 2) .* (param.c .* Z .* (((sign(P) + 1) ./ 2) .* P - param.P0 ) ./ (param.Hz + ((sign(P) + 1) ./ 2) .* P - param.P0))...
                       - P .* ((sign(P) + 1) ./ 2) .* param.r;

    dZdt = param.eps .* (((sign(Z) + 1) ./ 2) .* param.c .* Z .*( ((sign(P) + 1) ./ 2) .* P - param.P0 ) ./ (param.Hz + ((sign(P) + 1) ./ 2) .* P - param.P0))...
           - ((sign(Z) + 1) ./ 2) .*  param.d .* Z;

end


dydt = [dNdt; dPdt; dZdt];

end

