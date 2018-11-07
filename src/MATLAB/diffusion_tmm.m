% TMM
%
% Christian Ki??r and Anton Almgren
clear all
%close all
clc

%%
%load("../../bin/MITgcm/Matrix10/TMs/matrix_nocorrection.mat");
load('../../bin/MITgcm/Matrix5/TMs/matrix_nocorrection_annualmean.mat');
load('../../bin/MITgcm/grid.mat');
load('../../bin/MITgcm/config_data.mat');
load('../../bin/MITgcm/Matrix5/Data/boxes.mat')

%%

Ix = speye(nb,nb);
Aimpms = Aimpms^(72);%(24*60*60);

%finding where there is more than one layer
[fluxind1(:,1),fluxind1(:,2)]= find(bathy(:,:,1)==1);
[fluxind2(:,1),fluxind2(:,2)]= find(bathy(:,:,2)==1);
[val,ind1,ind2] = intersect(fluxind1,fluxind2,'rows','stable');

%fetching the diffusion values from Aimpms
diff_d = Aimpms(4448+ind2,ind1);
diff_d = diag(diff_d);%./(70-50)^2;
diff_d= full(diff_d);%*60*60*24/1200);

Dd = zeros(128,64);
Dd(:,:) = NaN;
for i = 1:length(diff_d)
    Dd(val(i,1),val(i,2)) = diff_d(i);
end

diff_u = Aimpms(ind1,4448+ind2);
diff_u = diag(diff_u);%./20^2;
diff_u= full(diff_u);%*60*60*24/1200);

Du = zeros(128,64);
Du(:,:) = NaN;
for i = 1:length(diff_u)
    Du(val(i,1),val(i,2)) = diff_u(i);
end
%%
Dd(isnan(Dd))=0;
figure
axesm eckert4; 
hold on
worldmap([y(1),y(end)],[x(1),x(end)])%,[y(1) y(end)],[x(1) x(end)])
% load coastlines
% plotm(coastlat,coastlon)
surfacem(y,x,Dd');
geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
colorbar
%caxis([0 0.6])
title('1->2')
%%
Du(isnan(Du))=0;
figure
axesm eckert4; 
hold on
worldmap([y(1),y(end)],[x(1),x(end)])%,[y(1) y(end)],[x(1) x(end)])
% load coastlines
% plotm(coastlat,coastlon)
surfacem(y,x,Du');
geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
colorbar
%caxis([0 0.6])
title('2->1')


