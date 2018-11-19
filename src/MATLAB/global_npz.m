% Global NPZ model
%
% Christian Kiaer and Anton Almgren
clear all
close all
clc
%% Load TMM
%load('../../bin/MITgcm/Matrix5/TMs/matrix_nocorrection_annualmean.mat');
load('../../bin/MITgcm/Matrix5/TMs/matrix_nocorrection_01.mat');
load('../../bin/MITgcm/grid.mat');
load('../../bin/MITgcm/config_data.mat');
load('../../bin/MITgcm/Matrix5/Data/boxes.mat')

%% Param setup
param.Hp = 0.5;   % Half-saturation constant Phyto
param.Hz = 1.0;   % Half
param.r = 0.035;   % Mortality rate
param.M = 50;     % Depth
param.D = 0.005;  % Diffusion

param.d = 0.035;   % Grazers loss
param.eps = 0.5;  % Grazing efficiency
param.P0 = 0.1;   % Grazing Threshold
param.N0 = 10;    % Deep nutrients

param.c = 1.0;    % Maximum Grazing

param.phi = 47;   % Test phi

loadPath = '../../bin/MITgcm/Matrix5/TMs/matrix_nocorrection_0';
loadPath1 =  '../../bin/MITgcm/Matrix5/TMs/matrix_nocorrection_';

%% Initial conditions
N = zeros(128,64,15);
P = zeros(128,64,15);
Z = zeros(128,64,15);

N(:,:,:) = 10;
P(:,:,1) = 1;
Z(:,:,1) = 0.1;

N = gridToMatrix(N, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
P = gridToMatrix(P, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
Z = gridToMatrix(Z, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');

Ix = speye(nb,nb);
Aexp = Ix + (12*60*60)*Aexp;
Aimp = Aimp^(36);
%%
surf_ind = length(find(bathy(:,:,1) == 1));
%%
%opts = odeset('NonNegative',1);

month = 1;
for i=1:200
    if mod(i, 61) == 0
        if month < 10
            load(strcat(loadPath, num2str(month), '.mat'));
        else
            load(strcat(loadPath1, num2str(month), '.mat'));
        end
        Aexp = Ix + (12*60*60)*Aexp;
        Aimp = Aimp^(36);
        month = month + 1;
        if month >= 13;
            month = 1;
        end
    end
    N =  Aimp * ( Aexp  * N);
    P =  Aimp * ( Aexp  * P);
    Z =  Aimp * ( Aexp  * Z);
    [t, Y] = ode23(@ode_npz, [0 0.5], [N;P;Z], [], param, i, surf_ind, Ybox(1:surf_ind));
    N = Y(end, 1:52749)'; 
    P = Y(end, 52749+1:52749*2)';
    Z = Y(end, 52749*2+1:52749*3)';
    N(N < 0) = 0;
    P(P < 0) = 0;
    Z(Z < 0) = 0;
    if mod(i,100) == 0
        i
%         Nn = matrixToGrid(N, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
%         Pn = matrixToGrid(P, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
%         Zn = matrixToGrid(Z, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
%         figure
%         surface(x,y, Nn(:,:,1)')
%         colorbar
%         %caxis([0,4])
%         figure
%         surface(x,y, Pn(:,:,1)')
%         colorbar
%         figure
%         surface(x,y, Zn(:,:,1)')
%         colorbar
    end
end
%%
N = matrixToGrid(N, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
P = matrixToGrid(P, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
Z = matrixToGrid(Z, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
%%
g(1:surf_ind,1:365) = exp(-(0.025) * param.M).*(1-0.8*sin(pi.*Ybox(1:surf_ind)/180).*cos(2.*pi.*(1:365)./365));
surface(g')
shading flat
%%
figure
surface(x,y, N(:,:,1)')
colorbar
%caxis([0,4])
figure
surface(x,y, P(:,:,1)')
colorbar
%caxis([0, 1])
figure
surface(x,y, Z(:,:,1)')
colorbar

