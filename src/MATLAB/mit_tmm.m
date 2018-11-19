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

%%1
surface(x,y,bathy(:,:,1)')
shading flat
C = zeros(128,64,15);

bathy(60,30,1)
C(:,:,1) = 1;

Ib=find(izBox>0);
nbb=length(Ib);
I = linspace(1,52749, 52749)';
Ia = linspace(1, 128*64*15, 128*64*15)';
mat = gridToMatrix(C, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
C = matrixToGrid(mat, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
figure
contourf(C(:,:,1)', 5)
colorbar
caxis([0 10])


%%
Ix = speye(nb,nb);
Aexpms = Ix + (12*60*60)*Aexpms;
Aimpms = Aimpms^(36);
%%
H = 2; %years
D = (1/2)/(730*H); 
Cn = mat;

for i=1:730*10%365*6
    Cn = matrixToGrid(Cn, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
    Cn(:,:,1) = 1;
    Cn = gridToMatrix(Cn, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');
    Cn =  Aimpms * ( Aexpms  * Cn) - Cn * D;

end
%Cn = Cn ;
%% Conservation of mass
a = volb/sum(volb);
for i = 1:1000:52000
    COM_c(i:i+1000) = a(i:i+1000).*mat(i:i+1000);
    COM_cn(i:i+1000) = a(i:i+1000).*Cn(i:i+1000);

end
COM_c(52000:52749)=a(52000:52749).*mat(52000:52749);
COM_cn(52000:52749)=a(52000:52749).*Cn(52000:52749);

COM = sum(COM_c)-sum(COM_cn);

Cn = matrixToGrid(Cn, [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat');


%%

figure
contourf(Cn(55:65,20:40,1)', 3)
colorbar
contourf(Cn(:,:,1)', 5)
colorbar
%caxis([0 100])
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

