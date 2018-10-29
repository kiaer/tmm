% TMM
%
% Christian Kiær and Anton Almgren
clear all
close all
clc

%%
load("../../bin/MITgcm/Matrix5/TMs/matrix_nocorrection_annualmean.mat");
load("../../bin/MITgcm/grid.mat");
loqad("../../bin/MITgcm/config_data.mat");
%%
surface(x,y,bathy(:,:,1)')
shading flat
C = zeros(128,64,15);

bathy(60,30,1)
C(60,30,1) = 1;

I = linspace(1,52749, 52749)';

mat = gridToMatrix(C, I, "../../bin/MITgcm/Matrix5/Data/boxes.mat", "../../bin/MITgcm/grid.mat");
%%
Cn = Aexpms* mat + Aimpms * mat;
