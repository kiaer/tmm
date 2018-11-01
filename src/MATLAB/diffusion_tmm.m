% TMM
%
% Christian Ki??r and Anton Almgren
clear all
close all
clc

%%
%load("../../bin/MITgcm/Matrix10/TMs/matrix_nocorrection.mat");
load('../../bin/MITgcm/Matrix5/TMs/matrix_nocorrection_annualmean.mat');
load('../../bin/MITgcm/grid.mat');
load('../../bin/MITgcm/config_data.mat');
load('../../bin/MITgcm/Matrix5/Data/boxes.mat')

%%

%finding where there is more than one layer
[fluxind1(:,1),fluxind1(:,2)]= find(bathy(:,:,1)==1);
[fluxind2(:,1),fluxind2(:,2)]= find(bathy(:,:,2)==1);
[val,ind1,ind2] = intersect(fluxind1,fluxind2,'rows','stable');

%fetching the diffusion values from Aimpms
diff_d = Aimpms(4448+ind2,ind1);
diff_d = diag(diff_d);
diff_d= full(diff_d);

Dd = zeros(128,64);
Dd(:,:) = NaN;
for i = 1:length(diff_d)
    Dd(val(i,1),val(i,2)) = diff_d(i);
end

diff_u = Aimpms(ind1,4448+ind2);
diff_u = diag(diff_u);
diff_u= full(diff_u);

Du = zeros(128,64);
Du(:,:) = NaN;
for i = 1:length(diff_u)
    Du(val(i,1),val(i,2)) = diff_u(i);
end
%%

figure
worldmap([y(1) y(end)],[x(1) x(end)])
hold on
% load coastlines
% plotm(coastlat,coastlon)
geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
surfacem(y,x,Dd');
title('down')
%%
figure
surface(x,y,Du');
title('up')





