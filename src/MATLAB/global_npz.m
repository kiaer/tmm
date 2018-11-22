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
param.r = 0.07;   % Mortality rate
param.M = 50;     % Depth
param.D = 0.005;  % Diffusion

param.d = 0.07;   % Grazers loss
param.eps = 0.5;  % Grazing efficiency
param.P0 = 0.1;   % Grazing Threshold
param.N0 = 10;    % Deep nutrients

param.c = 1.0;    % Maximum Grazing

param.phi = 47;   % Test phi

loadPath = '../../bin/MITgcm/Matrix5/TMs/matrix_nocorrection_0';
loadPath1 =  '../../bin/MITgcm/Matrix5/TMs/matrix_nocorrection_';

%% Initial conditions
% N = zeros(128,64,15);
% P = zeros(128,64,15);
% Z = zeros(128,64,15);
% 
% N(:,:,:) = 10;
% P(:,:,:) = 1;
% Z(:,:,:) = 0.1;

load('../../bin/init_values.mat')

N = gridToMatrix(N, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
P = gridToMatrix(P, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
Z = gridToMatrix(Z, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');

Ix = speye(nb,nb);
Aexp = Ix + (12*60*60)*Aexp;
Aimp = Aimp^(36);
%%
surf_ind = length(find(bathy(:,:,1) == 1));
%%
opts = odeset('RelTol',1e-2, 'AbsTol', 1e-5);%,'Stats','on');
opts = [];
month = 1;
layer = permute(sum(sum(bathy)),[3,1,2]);
layerind = [0 ; cumsum(layer) ];

% %% half layers
% layer1 = permute(sum(sum(bathy)),[3,1,2]);

% layer = zeros(2*length(layer1),1);
% layer(1:2:end-1) = ceil(layer1/2);
% layer(2:2:end) = floor(layer1/2);
% layerind = [0 ; cumsum(layer) ];
%%
for i=1:730*3
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

%     for j=1:length(N)
%         [t, Y] = ode23tb(@ode_npz, [0 0.5], [N(j) P(j) Z(j)], opts, param, i, surf_ind, Ybox(1:surf_ind), j);
%         if mod(j,10000) == 0
%             j
%         end
%         N(j) = Y(1);
%         P(j) = Y(2);
%         Z(j) = Y(3);
%     end
    for j=1:length(layer)
        [t, Y] = ode23(@ode_npz, [0 0.5], [N(layerind(j)+1:layerind(j+1))'...
             P(layerind(j)+1:layerind(j+1))' Z(layerind(j)+1:layerind(j+1))']...
            , opts, param, i, surf_ind, Ybox(1:surf_ind), j, layerind, layer);
        %if mod(j,1) == 0
        %    j
        %end
        N(layerind(j)+1:layerind(j+1)) = Y(end,1:layer(j));
        P(layerind(j)+1:layerind(j+1)) = Y(end,layer(j)+1:2*layer(j));
        Z(layerind(j)+1:layerind(j+1)) = Y(end,2*layer(j)+1:3*layer(j));
    end
    %N = Y(end, 1:52749)'; 
    %P = Y(end, 52749+1:52749*2)';
    %Z = Y(end, 52749*2+1:52749*3)';

    N(N < 0) = 0;
    P(P < 0) = 0;
    Z(Z < 0) = 0;
    if mod(i,60) == 0
        i
        Nn = matrixToGrid(N, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
        Pn = matrixToGrid(P, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
        Zn = matrixToGrid(Z, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
        figure(1)
        clf(1)
        surface(x,y, Nn(:,:,1)')
        %set(gcf, "zdata", Nn(:,:,1)')
        colorbar
        drawnow
        %caxis([0,4])
        figure(2)
        clf(2)
        surface(x,y, Pn(:,:,1)')
        colorbar
        drawnow
        figure(3)
        clf(3)
        surface(x,y, Zn(:,:,1)')
        colorbar
        drawnow
    end
end
%%
N = matrixToGrid(N, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
P = matrixToGrid(P, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
Z = matrixToGrid(Z, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
%%
save('../../bin/init_values.mat', 'N', 'P', 'Z')
%%
g(1:surf_ind,1:365) = exp(-(0.025) * param.M).*(1-0.8*sin(pi.*Ybox(1:surf_ind)/180).*cos(2.*pi.*(1:365)./365));
mon = [0 31 28 31 30 31 30 31 31 30 31 30 ];
ticks = 1+cumsum(mon);
labels = [{'J'} {'F'} {'M'} {'A'} {'M'} {'J'} {'J'} {'A'} {'S'} {'O'} {'N'} {'D'}];

surface(1:365,Ybox(1:surf_ind),g)
shading flat
xticks(ticks)
xticklabels(labels)
ylabel('Latitude')
c=colorbar;
c.Label.String  = 'growthrate';
xlim([1 sum(mon)+31])
ylim([-80 80])

%%
figure
surface(x,y, N(:,:,1)')
colorbar
title('Nutrients')
%caxis([0,4])
figure
surface(x,y, P(:,:,1)')
colorbar
%caxis([0, 1])  
figure
surface(x,y, Z(:,:,1)')
colorbar

%%
%N(isnan(N))=0;
Np = [N(:,:,:); N(1,:,:)];
xp = [x-x(1) ;360];

%%
figure
subplot(3,1,1)
hold on
axesm eckert4;
ax = worldmap('world');
setm(ax, 'Origin', [0 200 0])
surfacem(y,xp,N(:,:,1)');
geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5],'EdgeColor',[0.5 1.0 0.5]);
c=colorbar;z
c.Label.String='concentration';
title('Nutrients, Jan')

subplot(3,1,2)
hold on
axesm eckert4;
ax = worldmap('world');
setm(ax, 'Origin', [0 200 0])
surfacem(y,xp,P(:,:,1)');
geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5],'EdgeColor',[0.5 1.0 0.5]);
c=colorbar;
c.Label.String='concentration';
title('Phytoplankton, Jan')

subplot(3,1,3)
hold on
axesm eckert4;
ax = worldmap('world');
setm(ax, 'Origin', [0 200 0])
surfacem(y,xp,Z(:,:,1)');
geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5],'EdgeColor',[0.5 1.0 0.5]);
c=colorbar;
c.Label.String='concentration';
title('Zooplankton, Jan')

