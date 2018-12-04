% TMM
%
% Christian Ki??r and Anton Almgren
clear all
close all
clc

%%
%load("../../bin/MITgcm/Matwwwwrix10/TMs/matrix_nocorrection.mat");

load('../../bin/MITgcm/grid.mat');
load('../../bin/MITgcm/config_data.mat');
load('../../bin/MITgcm/Matrix5/Data/boxes.mat')

%%

%finding where there is more than one layer
[fluxind1(:,1),fluxind1(:,2)]= find(bathy(:,:,1)==1);
[fluxind2(:,1),fluxind2(:,2)]= find(bathy(:,:,2)==1);
[val,ind1,ind2] = intersect(fluxind1,fluxind2,'rows','stable');

num = [{'01'},{'02'},{'03'},{'04'},{'05'},{'06'},{'07'},{'08'},{'09'},...
    {'10'},{'11'},{'12'}];


for i=1:12
    if i>=10
        load(['../../bin/MITgcm/Matrix5/TMs/matrix_nocorrection_', num2str(i) ,'.mat']);
    else
        load(['../../bin/MITgcm/Matrix5/TMs/matrix_nocorrection_0', num2str(i) ,'.mat']);
    end
        

    Aimp = Aimp^(72);% diffusion [1/d]

    %fetching the diffusion values from Aimpms
    diff_d(:,i) = full(diag(Aimp(4448+ind2,ind1)));
   

    Dd(:,:,i) = zeros(128,64);
    Dd(:,:,i) = NaN;
    for j = 1:length(diff_d(:,i))
        Dd(val(j,1),val(j,2),i) = diff_d(j,i);
    end

    diff_u(:,i) = full(diag(Aimp(ind1,4448+ind2)));


    Du(:,:,i) = zeros(128,64);
    Du(:,:,i) = NaN;
    for j = 1:length(diff_u(:,i))
        Du(val(j,1),val(j,2)) = diff_u(j,i);
    end
end
%%
Dd(isnan(Dd))=0;
Ddp = [Dd(:,:,:); Dd(1,:,:)];
xp = [x-x(1) ;360];
%%
%xp = xp - 90
%xp = circshift(xp, -180)
mon ={ 'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
%setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);
%%
figure('Position', [0, 0, 700, 400]);
set(gcf,'color','w');
%subplot(2,2,1)
hold on
ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', 'Grid', 'on', 'Frame', 'on',...
    'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2)
ax.XColor = 'white'
ax.YColor = 'white'

%plabel('PlabelLocation',20, 'PLabelMeridian', 91)
surfacem(y,xp ,Ddp(:,:,1)');
%shading interp
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
caxis([0 0.6])
c = colorbar('southoutside', 'FontSize',14);
c.Label.String  = 'Diffusion rate [1/d]';
box off
print('../../fig/diff_jan', '-dpng', '-r300');

figure('Position', [0, 0, 700, 400]);
set(gcf,'color','w');
%subplot(2,2,1)
hold on
ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', 'Grid', 'on', 'Frame', 'on',...
    'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2)
ax.XColor = 'white'
ax.YColor = 'white'

%plabel('PlabelLocation',20, 'PLabelMeridian', 91)
surfacem(y,xp ,Ddp(:,:,4)');
%shading interp
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
caxis([0 0.6])
c = colorbar('southoutside', 'FontSize',14);
c.Label.String  = 'Diffusion rate [1/d]';
box off
pause(1)
print('../../fig/diff_apr', '-dpng', '-r300');


figure('Position', [0, 0, 700, 400]);
set(gcf,'color','w');

%subplot(2,2,1)
hold on
ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', 'Grid', 'on', 'Frame', 'on',...
    'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2)
ax.XColor = 'white'
ax.YColor = 'white'

%plabel('PlabelLocation',20, 'PLabelMeridian', 91)
surfacem(y,xp ,Ddp(:,:,7)');
%shading interp
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
caxis([0 0.6])
c = colorbar('southoutside', 'FontSize',14);
c.Label.String  = 'Diffusion rate [1/d]';
box off
pause(1)
print('../../fig/diff_jul', '-dpng', '-r300');

figure('Position', [0, 0, 700, 400]);
set(gcf,'color','w');

%subplot(2,2,1)
hold on
ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', 'Grid', 'on', 'Frame', 'on',...
    'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2)
ax.XColor = 'white'
ax.YColor = 'white'

%plabel('PlabelLocation',20, 'PLabelMeridian', 91)
surfacem(y,xp ,Ddp(:,:,10)');
%shading interp
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
caxis([0 0.6])
c = colorbar('southoutside', 'FontSize',14);
c.Label.String  = 'Diffusion rate [1/d]';
box off
pause(1)
print('../../fig/diff_oct', '-dpng', '-r300');

%%
% figure
% for i = 7%7:12
%     %subplot(2,3,i-6)
%     hold on
%     axesm eckert4;
%     h = worldmap('world');
%     set(h, 'Visible', 'off')
%     % load coastlines
%     % plotm(coastlat,coastlon)
%     surfacem(y,x,Dd(:,:,i)');
%     geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5],'EdgeColor',[0.5 1.0 0.5]);
%     colorbar
%     caxis([0 0.6])
%     title(mon(i))
% %     
%  end
%%
h = figure('Position', [50, 50, 900, 550]);
set(gcf,'color','w');

for i = 1:12
    
    %subplot(2,2,1)
    
    ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', 'Grid', 'on', 'Frame', 'on',...
        'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2);
    ax.XColor = 'white';
    ax.YColor = 'white';
    axis tight manual

    %plabel('PlabelLocation',20, 'PLabelMeridian', 91)
    surfacem(y,xp ,Ddp(:,:,i)');
    %shading interp
    geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
    caxis([0 0.6])
    c = colorbar('southoutside', 'FontSize',14);
    c.Label.String  = 'Diffusion rate [1/d]';
    box off
    title(mon(i))
    drawnow
    
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,'../../fig/downwelling.gif','gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,'../../fig/downwelling.gif','gif','WriteMode','append');
    end
end
%%
% Du(isnan(Du))=0;
% figure
% axesm eckert4; 
% hold on
% worldmap([y(1),y(end)],[x(1),x(end)])%,[y(1) y(end)],[x(1) x(end)])
% % load coastlines
% % plotm(coastlat,coastlon)
% surfacem(y,x,Du');
% geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
% colorbar
% %caxis([0 0.6])
% title('2->1')

%% original gif
% for i = 1:12
%     
%     h=figure;
%     hold on
%     axesm eckert4;
%     ax = worldmap('world');
%     setm(ax, 'Origin', [0 200 0])
%     surfacem(y,x,Dd(:,:,i)');
%     geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5],'EdgeColor',[0.5 1.0 0.5]);
%     caxis([0 0.6])
%     c=colorbar;
%     c.Label.String='[d^{-1}]';
%     title(mon(i))
%     
%     
%     % Capture the plot as an image
%     frame = getframe(h);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     % Write to the GIF File
%     if i == 1
%         imwrite(imind,cm,'downwelling','gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,'downwelling','gif','WriteMode','append');
%     end
%     
% end

