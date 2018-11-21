% TMM
%
% Christian Ki??r and Anton Almgren
clear all
close all
clc

%%
%load("../../bin/MITgcm/Matrix10/TMs/matrix_nocorrection.mat");

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
mon ={ 'Jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'};
    figure 
% for i = 1%:6
%     %subplot(2,3,i)
%     hold on
%     axesm eckert4;
%     ax = worldmap('world');
%     setm(ax, 'Origin', [0 200 0])
%     surfacem(y,xp,Dd(:,:,i)');
%     geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
%     caxis([0 0.6])
%     colorbar
%     title(mon(i))
%    
% end
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
for i = 1:12
    
    h=figure;
    hold on
    axesm eckert4;
    ax = worldmap('world');
    setm(ax, 'Origin', [0 200 0])
    surfacem(y,xp,Ddp(:,:,i)');
    geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5],'EdgeColor',[0.5 1.0 0.5]);
    caxis([0 0.6])
    c=colorbar;
    c.Label.String='[d^{-1}]';
    title(mon(i))
    
    
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,'downwelling','gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,'downwelling','gif','WriteMode','append');
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

