clear all;close all;clc;
tic

mkdir('Figures');

%% Parameter Definition
num_of_presynaptic_neurons = 500;
num_of_postsynaptic_neurons = 100;
p = 0.01; % Connectivity probabillity

fr = 10; % Firing rate

E=-70e-3; % Membrane resting potential
Rm=10e6; % Membrane resistance
tau_m=10e-3; % Membrane time constant
Vth=-54e-3; % Threshold before spike occurs
Vreset=-80e-3; % Potential after spike
Vspike = 20e-3; % Height of spike at AP
E_K = -80e-3;       % reversal potential for spike-rate adaptation (SRA) and 
                    % refractory current
                    
A = 4e-7; % This is the maximal size of the synaptic current.

pAC = 0.4;  % Probability of Associational Commisural connections
IACsize = 1e-3; % This is the current size of a single CA3 neuron. 

% Model parameters definition
U = 1e-4;
tauD = 120e-3;
tauF = 300e-3;


SW = 0.01; % This is the window size for the synchronicity window in msec
                    
                    
dt=0.0002;      
tmax=1;

%% Simulation
neuronGroup = simulate(num_of_presynaptic_neurons,num_of_postsynaptic_neurons,...
    p,fr,E,Rm,tau_m,Vth,Vreset,Vspike,E_K,dt,tmax,A,U,tauD,tauF, pAC, IACsize);

%% Plotting figures
V = neuronGroup.Post_Synaptic_Potential - E;
V(V == Vspike-E) = nan;

h = figure('Visible','Off');
hold all;
binWidth = 20;
bin = 0:binWidth:tmax*1e3;
val = nan(1,length(bin));
for i = 1 : size(V,1)
    [pks,locs] = findpeaks(V(i,:)*1e3); % This will locate the location and size of the EPSPs.
    if ~isempty(pks)
        locs = locs * dt * 1e3;
        locs(pks<=0) = [];
        pks(pks<=0) = [];
        
        plot(locs,pks,'.','MarkerSize',15);
        
%         locs(pks== 1e3*(Vspike - E)) = [];
%         pks(pks == 1e3*(Vspike - E)) =[];
        tmpks = pks;
        d = diff(locs); % finding the time interval between peaks
        Y = discretize(d,bin); % splitting the time intervals into the appropriate bins
        pks(1)= [];
        occ = mode(Y);
        curVal = nan(sum(Y==occ),length(bin));
        for j = 1 : nanmax(Y)
            if sum(Y==j) ~=0
                curVal(1:sum(Y==j),j) = pks(Y==j);
            end
        end
        val = [val ;curVal];
        pks = tmpks;
        
    end
    
end
xlabel('Time [msec]');
ylabel('EPSP [mV]');
title('EPSP as a function of time');

saveas(h,'Figures\EPSP as a function of time.bmp');

h = figure('Visible','Off');
hold all;
X = (0:length(nanmean(val))-1)*binWidth;
Y = nanmean(val);
Er = nanstd(val)/(sqrt(sum(~isnan(val(:)))));
bar(X,Y)
errorbar(X,Y,Er,'.');

% X = (0:length(nansum(val))-1)*binWidth;
% Y = nansum(val);
% bar(X,Y);

axis tight
xlabel('PSP interval width [msec]');
ylabel('PSP [mV]');
title('Average PSP as a function of interval width');
saveas(h,'Figures\PSP vs interval.bmp');


h = figure('Visible','Off');
V = neuronGroup.Post_Synaptic_Potential - E;
V(V<Vspike-E) = 0;
V(V == Vspike-E) = 1;
imagesc(~V)
colormap(gray);
xlabel('Time');
ylabel('Neuron');
title('Post Synaptic Raster Plot');
saveas(h,'Figures\Post Synaptic Raster Plot.bmp');

SyncIdx = nan(num_of_postsynaptic_neurons,1);
for i = 1 : num_of_postsynaptic_neurons
    idx = find(V(i,:));
    d = diff(idx)*dt;
    if ~isempty(sum(d<=SW))
        SyncIdx(i) = sum(d<=SW);
    end
end
SyncIdx(SyncIdx == 0) =nan;
disp(['The synchronicity index is: ' num2str(nanmean(SyncIdx))]);
    
    
    

toc