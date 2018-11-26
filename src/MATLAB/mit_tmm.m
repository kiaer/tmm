%% TMM MIT 2.8 degrees GCM
% This script is implemented to use the transport matrix method for
% transport of passive tracers in a global ocean.
% The method is described in detail here: https://github.com/samarkhatiwala/tmm
% Transport matrices and related data are available at 
% http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/
% 
% Authors: Christian Ki?r and Anton Almgren
% Fall 2018

clear all
close all
clc

%% Load transport matrices and configuration data

load('../../bin/MITgcm/Matrix5/TMs/matrix_nocorrection_annualmean.mat');
load('../../bin/MITgcm/grid.mat');
load('../../bin/MITgcm/config_data.mat');
load('../../bin/MITgcm/Matrix5/Data/boxes.mat')

%% Initializing arrays

C = zeros(128,64,15);                   % Concentration matrix for tracer
C(:,:,1) = 1;                           % Initial condition for tracer

% vector representation of C:
mat = gridToMatrix(C, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');


%% Preparing transport matrices for time stepping

Ix = speye(nb,nb);                      % Sparse Identity matrix 
Aexpms = Ix + (12*60*60)*Aexpms;        % Discretizing Aexpms for Explicit
                                        % euler time stepping using 12h
                                        % timesteps
Aimpms = Aimpms^(36);                   % modifying Aimpms to use 12h timesteps
%% C14 model
H = 5000;                               % years
D = (1/2)/(730*H);                      % decay
Cn = mat;                               % initializing concentration

% iterating through time for 10 yrs using 12h steps
for i=1:730*10
    % converting to grid format:
    Cn = matrixToGrid(Cn, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
    % setting surface to fixed value
    Cn(:,:,1) = 1;
    % Converting concentration field to matrix format (vector)
    Cn = gridToMatrix(Cn, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
    % time stepping with decay
    Cn =  Aimpms * ( Aexpms  * Cn) - Cn * D;
end
%% Conservation of mass
% only consistent with closed boundary conditions and no source terms
a = volb/sum(volb);         % relative volume
% calculating in small bits for matlab's sake
for i = 1:1000:52000
    COM_c(i:i+1000) = a(i:i+1000).*mat(i:i+1000);   % mass in each cell initially
    COM_cn(i:i+1000) = a(i:i+1000).*Cn(i:i+1000);   % mass in each cell after iteration
end
COM_c(52000:52749)=a(52000:52749).*mat(52000:52749);
COM_cn(52000:52749)=a(52000:52749).*Cn(52000:52749);

COM = sum(COM_c)-sum(COM_cn);       % Conservation of mass. Should be small.
%% Plotting

% converting concentration vector to grid
Cn = matrixToGrid(Cn, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');

% contour plot of surface values
figure
contourf(Cn(:,:,1)', 5)
colorbar

% Transect plot at longitude 30 (Atlantic)
figure
Cx = permute(Cn, [1,3,2]);
contourf(x,z,Cx(:,:,30)',5)
axis ij

%% 
figure
Cy = permute(Cn, [2,3,1]);
contourf(y(32:55),z,Cy(32:55,:,113)',6)
colormap(parula)
shading interp
colorbar
axis ij

