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

C = zeros(128,64,15);               % Concentration matrix for tracer
%C(60:80,12,1) = 10;                 % Initial condition for tracer Southern ocean
C(108:118,53) = 10;                 % Initial condition for tracer North atlantic

xp = [x-x(1) ;360];                 %x-vector for plotting

% vector representation of C:
mat = gridToMatrix(C, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
C = matrixToGrid(mat, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
%%
figure('Position', [0, 0, 700, 400]);
set(gcf,'color','w');
%subplot(2,2,1)
hold on
ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', 'Grid', 'on', 'Frame', 'on',...
    'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2);
ax.XColor = 'white';
ax.YColor = 'white';
box off

%plabel('PlabelLocation',20, 'PLabelMeridian', 91)
h=surfacem(y,xp ,C(:,:,1)');
%shading interp
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
c = colorbar('southoutside', 'FontSize',18);
c.Label.String  = 'Tracer Concentration';
set(gca,'FontSize', 14)
print('../../fig/tracer_atlantic_init', '-dpng', '-r300');



%% Preparing transport matrices for time stepping

Ix = speye(nb,nb);                      % Sparse Identity matrix 
Aexpms = Ix + (12*60*60)*Aexpms;        % Discretizing Aexpms for Explicit
                                        % euler time stepping using 12h
                                        % timesteps
Aimpms = Aimpms^(36);                   % modifying Aimpms to use 12h timesteps
%% Diagnostic plots
Cn = mat;                               % initializing concentration

% iterating through time for 5 yrs using 12h steps
for i=1:730*5
    Cn =  Aimpms * ( Aexpms  * Cn);
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


%%
% contour plot of surface values
figure('Position', [0, 0, 700, 400]);
set(gcf,'color','w');
%subplot(2,2,1)
hold on
ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', 'Grid', 'on', 'Frame', 'on',...
    'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2)
ax.XColor = 'white'
ax.YColor = 'white'

%plabel('PlabelLocation',20, 'PLabelMeridian', 91)
h=surfacem(y,xp ,Cn(:,:,1)');
%shading interp
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
c = colorbar('southoutside', 'FontSize',18);
c.Label.String  = 'Tracer Concentration';
%z1 = get(h,'ZData');
%set(h,'ZData',z1-10)  
%plot(-2:2, ones(5)*-1.08 ,'color','r','linewidth',1)

box off
print('../../fig/tracer_atlantic', '-dpng', '-r300');

%Transect plot at longitude 43w (Atlantic)
figure('Position', [0, 0, 700, 400]);
Cy = permute(Cn, [2,3,1]);
%xp1 = circshift(x, -90)
contourf(y(32:55),z, Cy(32:55,:,113)', 6)
xlabel('Latitude [\circ]')
ylabel('Depth [m]')
c = colorbar();
colormap(parula(7))
c.Label.String  = 'Tracer Concentration';
set(gca,'FontSize', 18)
set(gca,'Color',[0.8 0.8 0.8])

axis ij
print('../../fig/tracer_atlantic_transect', '-dpng', '-r300');


