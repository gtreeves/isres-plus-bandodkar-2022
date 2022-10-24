%sc_run
% This script runs the Dl/Cact model

clear
clc
close all


fhandle      = @calc_error_dc;
folderpath   = '/Volumes/Extreme/ISRES/Results/v3.3/dsat/';
filename     = 'Results_isres-plus_dsat_2022-06-08-21-04-05';



addpath('Plots')
filepath    = [folderpath,filename];
Filenames   = extractFileLocations(folderpath,'mat',true);
Filenames   = Filenames(contains(Filenames,filename));
total       = num2str(length(Filenames));

for i=1:length(Filenames)
    
    disp(['Running: ',num2str(i),'/',total]);
    matfile     = char(Filenames(i));    


    % load results and run model
    load(matfile,'xb','addParams','modelname');
    params               = [xb, addParams];
    [Fnew,penalty,~,Soln1] = fhandle(params,modelname);
    

    % plotting
    h = figure('position',[500,500,1450,500],'Visible','on'); 
    c = strsplit(matfile,{filesep,'.','_'}); 
    subplot(1,2,1)
    plotComparison(Soln1.dlNuc,Soln1.dlCactNuc,Soln1.T)    
    title(['Error = ',num2str(Fnew)],'fontsize',12)
    %makeTimeCourseMovie(Soln.totalDorsal)


end





