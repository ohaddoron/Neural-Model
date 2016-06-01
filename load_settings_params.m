function [settings, params] = load_settings_params(ModelMode)

%% Network_params
params.num_of_presynaptic_neurons = 3;
params.num_of_postsynaptic_neurons = 20;
params.p = 0.1; % Connectivity probabillity
params.fr = 25; % Average Firing rate

params.pAC = linspace(0,0.6,2);  % Probability of Associational Commisural connections
params.pAC = 0;
params.IACsize = 1e-10; % This is the current size of a single CA3 neuron. 

%% Pre-synaptic Activity params:

params.muISI = 25e-3;
params.muIBI = 500e-3;
params.sigmaISI = 3e-2;
params.sigmaIBI = 1.5 * params.sigmaISI;

%% neuronal params
params.refract = 10; % Refractory period length in ms
params.E=-56.9e-3; % Membrane resting potential
params.Rm=186e6; % Membrane resistance
params.tau_m=50e-3; % Membrane time constant
params.Vth=-42.3e-3; % Threshold before spike occurs
params.Vreset=-80e-3; % Potential after spike
params.Vspike = 20e-3; % Height of spike at AP
params.E_K = -80e-3;       % reversal potential for spike-rate adaptation (SRA) and 
                    % refractory current
%% synaptic parameters
params.A = 5e-8;  % This is the maximal size of the synaptic current of mossy-fibers.

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
settings.tmax=2;

settings.one_on_one=false; % a simple simulation of two neurons (pre and post)
settings.one_on_one_multi=false; 
settings.ISI_IBI = false;

if settings.one_on_one
    params.p = 1;
    params.num_of_presynaptic_neurons = 1;
    params.num_of_postsynaptic_neurons = 1;
end
settings.SW = 0.02; % This is the window size for the synchronicity window in msec
params.mu_lastAP = 200; % The mean value of the last firing time in the case of ISI and IBI firing if the neuron still hasn't fired.
params.sigma_lastAP = 100; % The std of the lasting firing time in the case of ISI and IBI firing if the neuron still hasn't fired.
settings.num_of_simulations = 10; % Number of simulations to run 
settings;
params;