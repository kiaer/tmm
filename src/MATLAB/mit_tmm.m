% TMM
%
% Christian Kiær and Anton Almgren
clear all
close all
clc

%%
%load("../../bin/MITgcm/Matrix10/TMs/matrix_nocorrection.mat");
load("../../bin/MITgcm/Matrix5/TMs/matrix_nocorrection_annualmean.mat");
load("../../bin/MITgcm/grid.mat");
load("../../bin/MITgcm/config_data.mat");
load('../../bin/MITgcm/Matrix5/Data/boxes.mat')
%%
surface(x,y,bathy(:,:,1)')
shading flat
C = zeros(128,64,15);

bathy(60,30,1)
C(60:90,10,1) = 10;
Ib=find(izBox>0);
nbb=length(Ib);
I = linspace(1,52749, 52749)';
Ia = linspace(1, 128*64*15, 128*64*15)';
mat = gridToMatrix(C, [], "../../bin/MITgcm/Matrix5/Data/boxes.mat", "../../bin/MITgcm/grid.mat");
C = matrixToGrid(mat, [], "../../bin/MITgcm/Matrix5/Data/boxes.mat", "../../bin/MITgcm/grid.mat");
figure
contourf(C(:,:,1)', 5)
colorbar
caxis([0 10])

%%
sum(mat)
Ix = speye(nb,nb);
Aexpms = Ix + 1200*Aexpms;
%%
Cn = mat;
for i=1:100000
    Cn =  Aimpms * ( Aexpms  * Cn);
    Cn( find(Cn < 0.000001)) = 0;
    %Cn = Cn / sum(Cn);
   % sum(Cn)
end
%Cn = Cn ;
sum(Cn, 'all', 'omitnan') 
Cn = matrixToGrid(Cn, [], "../../bin/MITgcm/Matrix5/Data/boxes.mat", "../../bin/MITgcm/grid.mat");
%%

figure
contourf(Cn(55:65,20:40,1)', 3)
colorbar
contourf(Cn(:,:,1)', 5)
colorbar
%caxis([0 100])
figure
Cx = permute(Cn, [1,3,2])
contourf(x,z,Cx(:,:,10)',5)
axis ij
