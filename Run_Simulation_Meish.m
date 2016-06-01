% There are two ways to run the simulation. 
% The first way is to use the varying synaptic parametrs which lead to
% strong depression, depression etc.
% The second way is to set a gradient of AC connection probabilities. 
% The method in which the simulation works is determined by the length of
% the pAC vector, meaning, if there is only 1 element in it, the simulation
% will use the first way, if it has more than 1 element, it will use the
% second way.
% The method in which the model runs will be determined using the data set
% in the load_settings_params function. In that function, the pAC vector is
% set 
%% Program start
clear all;close all;clc;
tic
ModelMode = 1;
[settings, params] = load_settings_params(ModelMode); % From this point, the simulation method is set.

%% Removing any previously existing folders
% In order for us not to get confused with any previously saved data, this
% part removes any old folders that existed in the current folder path.
try 
    rmdir(fullfile(pwd,'Figures'),'s');
catch
    
end

for i = 1 : 5
    try
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
        rmdir(path2,'s');
    catch
        continue;
    end
end

mkdir('Figures');
                                 
%% Simulation
num_of_simulations = settings.num_of_simulations;

for n = 1 : num_of_simulations
    
    
    neuronGroup = cell(numel(params.pAC),1);
    if numel(params.pAC) > 1 % If this is true, the simulation method is a gradient of AC connections
        ModelMode = 3;
        [settings, params] = load_settings_params(ModelMode);
        for i = 1 : numel(params.pAC)
            disp(strcat('Probability of AC connections: ',num2str(params.pAC(i))));
            neuronGroup{i} = simulate(params, settings,i);
            close all;
        end


    else % Otherwise, the varying synaptic parameters will run
        for i = 1 : 5
            ModelMode = i;
            [settings, params] = load_settings_params(ModelMode);
            switch i
                case 1
                    path2 = fullfile(pwd,'Strong Depression');
                    disp('Strong Depression mode');
                case 2
                    path2 = fullfile(pwd,'Depression');
                    disp('Depression mode');
                case 3
                    path2 = fullfile(pwd,'Facilitation-Depression');
                    disp('Facilitation-Depression mode');
                case 4
                    path2 = fullfile(pwd,'Facilitation');
                    disp('Facilitation mode');
                case 5
                    path2 = fullfile(pwd,'Strong Facilitation');
                    disp('Strong Facilitation mode');
            end
            outPath = path2;
            neuronGroup{i} = simulate(params,settings,1);
            copyfile(fullfile(pwd,'Figures'),path2);
            rmdir(fullfile(pwd,'Figures'),'s');
            mkdir(fullfile(pwd,'Figures'));
            close all;
        end
        neural_struct{n} = neuronGroup;
            
        
    
    end
    disp(['Simulation # ' num2str(n) ' Done!']);
end
   
%% Plotting figures + multi simulation plots
Sync_Indicator = [];
sub_treshold_EPSP = [];
num_of_post_synaptic_spikes = [];
for n = 1 : num_of_simulations   
    neuronGroup = neural_struct{n};
    [Sync_Indicator_tmp,num_of_post_synaptic_spikes_tmp,sub_treshold_EPSP_tmp] = plotFigures (neuronGroup,params,settings);
    Sync_Indicator = cat(2,Sync_Indicator,Sync_Indicator_tmp);
    num_of_post_synaptic_spikes= cat(2,num_of_post_synaptic_spikes,num_of_post_synaptic_spikes_tmp');
    sub_treshold_EPSP = cat(2,sub_treshold_EPSP,sub_treshold_EPSP_tmp');
end
    
    
h = figure('Visible','Off');
bar(nanmean(Sync_Indicator,2));
hold all;
errorbar(1:numel(nanmean(Sync_Indicator,2)),nanmean(Sync_Indicator,2),nanstd(Sync_Indicator,0,2));
if numel(params.pAC) > 1
    xlabel('Probability of AC connections');
    ylabel('Synchronicity indicator');
    title('Synchronicity indicator as a function of connection probability');
    set(gca,'XTickLabel',num2cell(params.pAC));
else
    xlabel('Model Mode');
    ModelModes = {'Strong Depression','Depression','Facilitation-Depression','Facilitation','Strong Facilitation'};
    set(gca,'XTickLabel',ModelModes,'XTickLabelRotation',45);
    title('Synchronicity indicator as a function of model mode');
    ylabel('Synchronicity indicator');
end
saveas(h,fullfile(pwd,'Figures','Multiple simulations sync indicator.bmp'));

h = figure('Visible','Off');
bar(nanmean(sub_treshold_EPSP,2))
hold all;
errorbar(1:numel(nanmean(sub_treshold_EPSP,2)),nanmean(sub_treshold_EPSP,2),nanstd(sub_treshold_EPSP,0,2));
axis tight;
if numel(params.pAC) > 1
    xlabel('Probability of AC connections');
    ylabel('Total amount of sub threshold EPSPs/Number of post synaptic neurons');
    title('Total amount of sub threshold EPSPs as a function of connection probability');
    set(gca,'XTickLabel',num2cell(params.pAC));
else
    xlabel('Model Mode');
    ModelModes = {'Strong Depression','Depression','Facilitation-Depression','Facilitation','Strong Facilitation'};
    set(gca,'XTickLabel',ModelModes,'XTickLabelRotation',45);
    title('Total amount of sub threshold EPSPs as a function of model mode');
    ylabel('Total amount of sub threshold EPSPs/Number of post synaptic neurons');

end
saveas(h,fullfile(pwd,'Figures','Multiple simulations amount of sub threshold EPSPs.bmp'));

h = figure('Visible','Off');
bar(nanmean(num_of_post_synaptic_spikes,2));
hold all;
errorbar(1:numel(nanmean(num_of_post_synaptic_spikes,2)),nanmean(num_of_post_synaptic_spikes,2),nanstd(num_of_post_synaptic_spikes,0,2));
if numel(params.pAC) > 1
    xlabel('Probability of AC connections');
    ylabel('Number of post synaptic spikes/Number of post synaptic neurons');
    title('Number of post synaptic spikes as a function of connection probability');
    set(gca,'XTickLabel',num2cell(params.pAC));
else
    xlabel('Model Mode');
    ModelModes = {'Strong Depression','Depression','Facilitation-Depression','Facilitation','Strong Facilitation'};
    set(gca,'XTickLabel',ModelModes,'XTickLabelRotation',45);
    title('Number of post synaptic spikes as a function of model mode');
    ylabel('Number of post synaptic spikes/Number of post synaptic neurons'); 
end
saveas(h,fullfile(pwd,'Figures','Multiple simulations number of post synaptic spikes.bmp'));



toc