clear all;close all;clc;
tic
ModelMode = 1;
[settings, params] = load_settings_params(ModelMode);

mkdir('Figures');
                                 
%% Simulation
neuronGroup = cell(numel(params.pAC),1);
if numel(params.pAC) > 1
    for i = 1 : numel(params.pAC)
        neuronGroup{i} = simulate(params, settings,i);
        close all;
    end
    plotFigures (neuronGroup,params,settings)

else
    for i = 1 : 5
        ModelMode = i;
        [settings, params] = load_settings_params(ModelMode);
        neuronGroup{i} = simulate(params,settings,1);
        close all;
        path1 = fullfile(pwd,'Figures');
        switch i
            case 1
                path2 = fullfile(pwd,'Strong Depression');
            case 2
                path2 = fullfile(pwd,'Depression');
            case 3
                path2 = fullfile(pwd,'Facilitation-Depression');
            case 4
                path2 = fullfile(pwd,'Facilitation');
            case 5
                path2 = fullfile(pwd,'Strong Facilitation');
        end
        plotFigures (neuronGroup(i),params,settings)
        copyfile(path1,path2);
        rmdir(fullfile(pwd,'Figures'),'s');
        mkdir(fullfile(pwd,'Figures'));
    end
    
end




toc