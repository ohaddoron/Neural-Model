function [settings, params] = load_settings_params(ModelMode)

%% Network_params
params.num_of_presynaptic_neurons = 100;
params.num_of_postsynaptic_neurons = 20;
params.p = 0.01; % Connectivity probabillity
params.fr = 1; % Average Firing rate

% params.pAC = linspace(0,0.2,5);  % Probability of Associational Commisural connections
params.pAC = 0;
params.IACsize = 1e-10; % This is the current size of a single CA3 neuron. 

%% Pre-synaptic Activity params:

params.muISI = 25e-3;
params.muB = 250e-3;
params.sigmaA = 3e-2;
params.sigmaB = 1.5 * params.sigmaA;

%% neuronal params
params.E=-70e-3; % Membrane resting potential
params.Rm=10e6; % Membrane resistance
params.tau_m=10e-3; % Membrane time constant
params.Vth=-54e-3; % Threshold before spike occurs
params.Vreset=-80e-3; % Potential after spike
params.Vspike = 20e-3; % Height of spike at AP
params.E_K = -80e-3;       % reversal potential for spike-rate adaptation (SRA) and 
                    % refractory current
%% synaptic parameters
params.A = 4e-7;  % This is the maximal size of the synaptic current of mossy-fibers.

switch ModelMode 
    case 1 % Strong depression
        params.tauD = 1.7;
        params.tauF = 0.02;
        params.U = 0.7;
        params.f = 0.05;
    case 2 % Depression
        params.tauD = 0.5;
        params.tauF = 0.05;
        params.U = 0.5;
        params.f = 0.05;
    case 3 % Facilitation-Depression
        params.tauD = 0.2;
        params.tauF = 0.2;
        params.U = 0.25;
        params.f = 0.3;
    case 4 % Facilitation
        params.tauD = 0.05;
        params.tauF = 0.5;
        params.U = 0.15;
        params.f = 0.15;
    case 5 % Strong Facilitation
        params.tauD = 0.02;
        params.tauF = 1.7;
        params.U = 0.1;
        params.f = 0.11;
end


%% settings
settings.dt=0.0002;      
settings.tmax=1;

settings.one_on_one=false; % a simple simulation of two neurons (pre and post)

if settings.one_on_one
    params.p = 1;
    params.num_of_presynaptic_neurons = 1;
    params.num_of_postsynaptic_neurons = 1;
end
settings.SW = 0.01; % This is the window size for the synchronicity window in msec


settings;
params;