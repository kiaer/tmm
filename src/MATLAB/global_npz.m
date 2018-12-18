%% Global NPZ model
% 
% Christian Kiaer and Anton Almgren
clear all
close all
clc
%% Load Initial January TM and configuations
load('../../bin/MITgcm/Matrix5/TMs/matrix_nocorrection_01.mat');
load('../../bin/MITgcm/grid.mat');
load('../../bin/MITgcm/config_data.mat');
load('../../bin/MITgcm/Matrix5/Data/boxes.mat')

% Load paths for switching TM
loadPath = '../../bin/MITgcm/Matrix5/TMs/matrix_nocorrection_0';
loadPath1 =  '../../bin/MITgcm/Matrix5/TMs/matrix_nocorrection_';

% Preparing for timestepping. Using 43200 sec timesteps (0.5 days)
Ix = speye(nb,nb);
Aexp = Ix + (12*60*60)*Aexp;
Aimp = Aimp^(36);

%% Param setup
param.Hp = 0.5;   % Half-saturation constant Phyto
param.Hz = 1.0;   % Half
param.r = 0.07;   % Mortality rate
param.M = 50;     % Depth

param.d = 0.07;   % Grazers loss
param.eps = 0.5;  % Grazing efficiency
param.P0 = 0.1;   % Grazing Threshold
param.N0 = 10;    % Deep nutrients

param.c = 1.0;    % Maximum Grazing
%% Initial conditions
% Uncomment if no spun up version is available
% N = zeros(128,64,15);
% P = zeros(128,64,15);
% Z = zeros(128,64,15);
% 
% N(:,:,:) = 10;
% P(:,:,:) = 1;
% Z(:,:,:) = 0.1;

% Loads initial values from previous runs. Comment if no data is available
load('../../bin/init_values.mat')

% Convert from grid form to matrix form
N = gridToMatrix(N, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
P = gridToMatrix(P, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
Z = gridToMatrix(Z, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');

month = 0;

% Toplayer monthly mean vector
Nm = zeros(4448,12);
Pm = zeros(4448,12);
Zm = zeros(4448,12);

% Toplayer half daily vector
Nd = zeros(4448,730);
Pd = zeros(4448,730);
Zd = zeros(4448,730);

mon = [0 31 28 31 30 31 30 31 31 30 31 30 ];
ticks = (1+cumsum(mon))*2;
ticks(1) = 1;
labels = [{'J'} {'F'} {'M'} {'A'} {'M'} {'J'} {'J'} {'A'} {'S'} {'O'} {'N'} {'D'}];
%% Misc
Np = [N(:,:,:); N(1,:,:)];
xp = [x-x(1) ;360];

% Finding the number of indicies in the surface layer
surf_ind = length(find(bathy(:,:,1) == 1));
    
% Split into layers
layer = permute(sum(sum(bathy)),[3,1,2]);
layerind = [0 ; cumsum(layer)];


%% Running the TM with the simple NPZ model
for i=1:730
    
    % Test for time to change monthly TM
    if ismember(i, ticks)
        month = month + 1;
        % Reset to Jan
        if month > 12;
            month = 1;
        end
        % Load TM
        if month < 10
            load(strcat(loadPath, num2str(month), '.mat'));
        else
            load(strcat(loadPath1, num2str(month), '.mat'));
        end
        
        % Preparing for timestepping. 43200s.
        Aexp = Ix + (12*60*60)*Aexp;
        Aimp = Aimp^(36);
     
    end
    N =  Aimp * ( Aexp  * N);
    P =  Aimp * ( Aexp  * P);
    Z =  Aimp * ( Aexp  * Z);

    for j=1:length(layer)
        [t, Y] = ode23(@ode_npz, [0 0.5], [N(layerind(j)+1:layerind(j+1))'...
             P(layerind(j)+1:layerind(j+1))' Z(layerind(j)+1:layerind(j+1))']...
            , [], param, i, surf_ind, Ybox(1:surf_ind), j, layer);
        
        N(layerind(j)+1:layerind(j+1)) = Y(end,1:layer(j));
        P(layerind(j)+1:layerind(j+1)) = Y(end,layer(j)+1:2*layer(j));
        Z(layerind(j)+1:layerind(j+1)) = Y(end,2*layer(j)+1:3*layer(j));
    end

    N(N < 0) = 0;
    P(P < 0) = 0;
    Z(Z < 0) = 0;
    
    Nm(:,month) = Nm(:,month) + N(1:surf_ind);
    Pm(:,month) = Pm(:,month) + P(1:surf_ind);
    Zm(:,month) = Zm(:,month) + Z(1:surf_ind);
    
    Nd(:,i) = N(1:surf_ind);
    Pd(:,i) = P(1:surf_ind);
    Zd(:,i) = Z(1:surf_ind);
    
    % 
    if mod(i,60) == 0
       disp(['t=',num2str(i)]);
    end
end
%% Convert back to grid format
N = matrixToGrid(N, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
P = matrixToGrid(P, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
Z = matrixToGrid(Z, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
%% Monthly mean plots
mon_mean = [31 28 31 30 31 30 31 31 30 31 30 31]*2;
labels = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};

%% ------------------------------------------------------------------------
%   Plot Production
%  ------------------------------------------------------------------------
%% NUTRIENTS
for i=1:3:12
    if i==10 
        h = figure('Position', [0, 0, 700, 400]);
    else
        h = figure('Position', [0, 0, 700, 350]);
    end
    set(gcf,'color','w');
    Nm1(:,:,i) = matrixToGrid(Nm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
%     Pm1(:,:,i) = matrixToGrid(Pm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
%     Zm1(:,:,i) = matrixToGrid(Zm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
    
     ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', 'Grid', 'on', 'Frame', 'on',...
        'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2);
    ax.XColor = 'white';
    ax.YColor = 'white';
    axis tight manual

    %plabel('PlabelLocation',20, 'PLabelMeridian', 91)
    surfacem(y,xp ,Nm1(:,:,i)');
    %shading interp
    geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
    if i==10
        caxis([0 10])
        c = colorbar('southoutside', 'FontSize',18);
        c.Label.String  = 'Concentration [mmol m^{-3}]';
    end
    box off
    %title('Nutrients')
    %text(0, 1, labels(i),'Units','normalized')
    
    print(['../../fig/global_nutrients_',num2str(i)], '-dpng', '-r300');
end
%% PHYTOPLANKTON
for i=1:3:12
    if i==10 
        h = figure('Position', [0, 0, 700, 400]);
    else
        h = figure('Position', [0, 0, 700, 350]);
    end
    set(gcf,'color','w');
%     Nm1(:,:,i) = matrixToGrid(Nm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
    Pm1(:,:,i) = matrixToGrid(Pm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
%     Zm1(:,:,i) = matrixToGrid(Zm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
    
     ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', 'Grid', 'on', 'Frame', 'on',...
        'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2);
    ax.XColor = 'white';
    ax.YColor = 'white';
    axis tight manual

    %plabel('PlabelLocation',20, 'PLabelMeridian', 91)
    surfacem(y,xp ,Pm1(:,:,i)');
    %shading interp
    geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
    if i ==10
        caxis([0 10])
        c = colorbar('southoutside', 'FontSize',18);
        c.Label.String  = 'Concentration [mmol m^{-3}]';
    end
    box off
    %title('Phytoplankton')
    
    print(['../../fig/global_phytoplankton_',num2str(i)], '-dpng', '-r300');
end

%% ZOOPLANKTON
for i=1:3:12
    h = figure('Position', [0, 0, 700, 400]);
    set(gcf,'color','w');
%     Nm1(:,:,i) = matrixToGrid(Nm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
%     Pm1(:,:,i) = matrixToGrid(Pm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
    Zm1(:,:,i) = matrixToGrid(Zm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
    
     ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', 'Grid', 'on', 'Frame', 'on',...
        'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2);
    ax.XColor = 'white';
    ax.YColor = 'white';
    axis tight manual

    %plabel('PlabelLocation',20, 'PLabelMeridian', 91)
    surfacem(y,xp ,Zm1(:,:,i)');
    %shading interp
    geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
    caxis([0 2])
    c = colorbar('southoutside', 'FontSize',18);
    c.Label.String  = 'Concentration [mmol m^{-3}]';
    box off
    %title('Zooplankton')
    %text(0, 1, labels(i),'Units','normalized')
    
    print(['../../fig/global_zooplankton_',num2str(i)], '-dpng', '-r300');
end
%% -------------------------------------------------------------------------
%   GIF PRODUCTION
%  -------------------------------------------------------------------------
% h = figure('Position', [50, 50, 900, 550]);
% set(gcf,'color','w');
% for i=1:12
% % h = figure('Position', [50, 50, 900, 550]);
% % set(gcf,'color','w');
%     Nm1(:,:,i) = matrixToGrid(Nm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
% %     Pm1(:,:,i) = matrixToGrid(Pm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
% %     Zm1(:,:,i) = matrixToGrid(Zm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
%     
%      ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', 'Grid', 'on', 'Frame', 'on',...
%         'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2);
%     ax.XColor = 'white';
%     ax.YColor = 'white';
%     axis tight manual
% 
%     %plabel('PlabelLocation',20, 'PLabelMeridian', 91)
%     surfacem(y,xp ,Nm1(:,:,i)');
%     %shading interp
%     geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
%     caxis([0 10])
%     c = colorbar('southoutside', 'FontSize',14);
%     c.Label.String  = 'Concentration [mmol m^{-3}]';
%     box off
%     title('Nutrients')
%     text(0, 1, labels(i),'Units','normalized')
%     
% %     print(['../../fig/nutrients_',num2str(i)], '-dpng', '-r300');
% 
%     drawnow
%     frame = getframe(h);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     % Write to the GIF File
%     if i == 1
%         imwrite(imind,cm,'../../fig/seasonal_nutrients.gif','gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,'../../fig/seasonal_nutrients.gif','gif','WriteMode','append');
%     end
%     clf
% end
% %% Phytoplankton GIF
% h = figure('Position', [50, 50, 900, 550]);
% set(gcf,'color','w');
% for i=1:12
% % h = figure('Position', [50, 50, 900, 550]);
% % set(gcf,'color','w');
%     %Nm1(:,:,i) = matrixToGrid(Nm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
%     Pm1(:,:,i) = matrixToGrid(Pm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
% %     Zm1(:,:,i) = matrixToGrid(Zm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
%     
%      ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', 'Grid', 'on', 'Frame', 'on',...
%         'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2);
%     ax.XColor = 'white';
%     ax.YColor = 'white';
%     axis tight manual
% 
%     %plabel('PlabelLocation',20, 'PLabelMeridian', 91)
%     surfacem(y,xp ,Pm1(:,:,i)');
%     %shading interp
%     geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
%     caxis([0 10])
%     c = colorbar('southoutside', 'FontSize',14);
%     c.Label.String  = 'Concentration [mmol m^{-3}]';
%     box off
%     title('Phytoplankton')
%     text(0, 1, labels(i),'Units','normalized')
% 
% %    print(['../../fig/seasonal_phytoplankton_',num2str(i)], '-dpng', '-r300');
% 
%     drawnow
%     frame = getframe(h);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     % Write to the GIF File
%     if i == 1
%         imwrite(imind,cm,'../../fig/seasonal_phytoplankton.gif','gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,'../../fig/seasonal_phytoplankton.gif','gif','WriteMode','append');
%     end
%     clf
% end
% %% Zooplankton GIF
% h = figure('Position', [50, 50, 900, 550]);
% set(gcf,'color','w');
% for i=1:12
% %h = figure('Position', [50, 50, 900, 550]);
% %set(gcf,'color','w');
%     %Nm1(:,:,i) = matrixToGrid(Nm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
% %     Pm1(:,:,i) = matrixToGrid(Pm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
%     Zm1(:,:,i) = matrixToGrid(Zm(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat') / mon_mean(i);
%     
%      ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', 'Grid', 'on', 'Frame', 'on',...
%         'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2);
%     ax.XColor = 'white';
%     ax.YColor = 'white';
%     axis tight manual
% 
%     %plabel('PlabelLocation',20, 'PLabelMeridian', 91)
%     surfacem(y,xp ,Zm1(:,:,i)');
%     %shading interp
%     geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
%     caxis([0 4])
%     c = colorbar('southoutside', 'FontSize',14);
%     c.Label.String  = 'Concentration [mmol m^{-3}]';
%     box off
%     title('Zooplankton')
%     text(0, 1, labels(i),'Units','normalized')
% 
% %    print(['../../fig/seasonal_zooplankton_',num2str(i)], '-dpng', '-r300');
%     
%     drawnow
%     frame = getframe(h);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     % Write to the GIF File
%     if i == 1
%         imwrite(imind,cm,'../../fig/seasonal_zoo.gif','gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,'../../fig/seasonal_zoo.gif','gif','WriteMode','append');
%     end
%     clf
% end

%% Plot transect of the Atlantic
for i=1:730
    Nd1(:,:,i) = matrixToGrid(Nd(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
    Pd1(:,:,i) = matrixToGrid(Pd(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
    Zd1(:,:,i) = matrixToGrid(Zd(:,i), (1:surf_ind), '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
end
%%
Pd1p = permute(Pd1, [2,3,1]);
Nd1p = permute(Nd1, [2,3,1]);
Zd1p = permute(Zd1, [2,3,1]);
days = 0.5:0.5:365;
figure
surface(days,y,Pd1p(:,:,121));
xlim([0 365])
ylim([-71 71])
shading interp
c = colorbar('FontSize',14)
c.Label.String  = 'Concentration [mmol m^{-3}]';
xlabel('Time [days]')
ylabel('Latitude [\circ]')
print('../../fig/atlantic_phyto','-dpng', '-r300')
%%
figure
surface(days,y,Zd1p(:,:,121));
shading interp
figure
surface(days,y,Nd1p(:,:,121));
shading interp
%% Save new initial values
%save("../../bin/init_values.mat", 'N', 'P', 'Z')
%%
g(1:surf_ind,1:365) = exp(-(0.025) * param.M).*(1-0.8*sin(pi.*Ybox(1:surf_ind)/180).*cos(2.*pi.*(1:365)./365));

figure
surface(1:365,Ybox(1:surf_ind),g)
shading interp
xticks(ticks /2 + 1)
xticklabels(labels)
ylabel('Latitude [\circ]')
xlabel('Month')
c=colorbar;
c.Label.String  = 'Growthrate';
xlim([1 sum(mon)+31])
ylim([-77.3438 80.1562])
print('../../fig/growthrate', '-dpng', '-r300');

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
h = figure('rend','painters','pos',[10 10 900 1200])
subplot(3,1,1)
hold on
axesm eckert4;
ax = worldmap('world');
setm(ax, 'Origin', [0 200 0])
surfacem(y,xp,N(:,:,1)');
geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5],'EdgeColor',[0.5 1.0 0.5]);
c=colorbar;
c.Label.String='concentration';
title(strcat('Nutrients month = ', num2str(month)))

subplot(3,1,2)
hold on
axesm eckert4;
ax = worldmap('world');
setm(ax, 'Origin', [0 200 0])
surfacem(y,xp,P(:,:,1)');
geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5],'EdgeColor',[0.5 1.0 0.5]);
c=colorbar;
c.Label.String='concentration';
title(strcat('Phytoplankton month = ', num2str(month)))

subplot(3,1,3)
hold on
axesm eckert4;
ax = worldmap('world');
setm(ax, 'Origin', [0 200 0])
surfacem(y,xp,Z(:,:,1)');
geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5],'EdgeColor',[0.5 1.0 0.5]);
c=colorbar;
c.Label.String='concentration';
title(strcat('Zooplankton month = ', num2str(month)))