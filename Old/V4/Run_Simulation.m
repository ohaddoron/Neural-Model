clear all;close all;clc;
tic

mkdir('Figures');

%% Parameter Definition
num_of_presynaptic_neurons = 500;
num_of_postsynaptic_neurons = 100;
p = 0.01; % Connectivity probabillity

fr = 1; % Firing rate

E=-70e-3; % Membrane resting potential
Rm=10e6; % Membrane resistance
tau_m=10e-3; % Membrane time constant
Vth=-54e-3; % Threshold before spike occurs
Vreset=-80e-3; % Potential after spike
Vspike = 20e-3; % Height of spike at AP
E_K = -80e-3;       % reversal potential for spike-rate adaptation (SRA) and 
                    % refractory current
                    
A = 1e-8; % This is the maximal size of the synaptic current.

pAC = linspace(0,0.4,5);  % Probability of Associational Commisural connections
% pAC = 0;
IACsize = 1e-3; % This is the current size of a single CA3 neuron. 

% Model parameters definition
% U = 1e-4;
% tauD = 120e-3;
% tauF = 300e-3;
ModelMode = 5;
switch ModelMode 
    case 1 % Strong depression
        tauD = 1.7;
        tauF = 0.02;
        U = 0.7;
        f = 0.05;
    case 2 % Depression
        tauD = 0.5;
        tauF = 0.05;
        U = 0.5;
        f = 0.05;
    case 3 % Facilitation-Depression
        tauD = 0.2;
        tauF = 0.2;
        U = 0.25;
        f = 0.3;
    case 4 % Facilitation
        tauD = 0.05;
        tauF = 0.5;
        U = 0.15;
        f = 0.15;
    case 5 % Strong Facilitation
        tauD = 0.02;
        tauF = 1.7;
        U = 0.1;
        f = 0.11;
end



SW = 0.01; % This is the window size for the synchronicity window in msec
                    
                    
dt=0.0002;      
tmax=1;

%% Simulation
neuronGroup = cell(numel(pAC),1);
for i = 1 : numel(pAC)
    neuronGroup{i} = simulate(num_of_presynaptic_neurons,num_of_postsynaptic_neurons,...
        p,fr,E,Rm,tau_m,Vth,Vreset,Vspike,E_K,dt,tmax,A,U,tauD,tauF, pAC(i), IACsize , f);
    close all;
end

plotFigures (neuronGroup,pAC,E,Vspike,tmax,dt,num_of_postsynaptic_neurons,SW)


toc