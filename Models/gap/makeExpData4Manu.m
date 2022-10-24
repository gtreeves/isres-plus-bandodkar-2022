clear
clc
close all

folderpath = 'Models/gap/';

load([folderpath,'Mats/gap_gene_nc13.mat'])
load([folderpath,'Mats/TC_data_3592.mat'])
load([folderpath,'Mats/gapgrn.mat'])


% Only need to get nc 13 data!
x13  = linspace(params.xL,params.xU,params.M13)';
Hb13 = interp1(exp13{1}.x,exp13{1}.y,x13,'linear','extrap');
Kr13 = interp1(exp13{2}.x,exp13{2}.y,x13,'linear','extrap');
Gt13 = interp1(exp13{3}.x,exp13{3}.y,x13,'linear','extrap');
NC13 = [Hb13, Kr13, Gt13];


NC14 = TCData;


save([folderpath,'Mats/ExpDataManu.mat'],'NC13','NC14')

