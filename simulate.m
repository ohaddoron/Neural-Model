function neuronGroup = simulate(params, settings,n)

% Leaky-integrate-and-fire neuron with spike rate adaptation and a
% refractory current.





Connectivity_Ratio = connect(params.num_of_presynaptic_neurons, params.num_of_postsynaptic_neurons, params.p);  % Drawing the connectivity map using the connectivity probability
neuronGroup.Connectivity_Ratio = Connectivity_Ratio; 
neuronGroup.num_of_presynaptic_neurons = params.num_of_presynaptic_neurons;
neuronGroup.num_of_postsynaptic_neurons = params.num_of_postsynaptic_neurons;

AC_Connectivity_Ratio = connect(params.num_of_postsynaptic_neurons, params.num_of_postsynaptic_neurons, params.pAC(n)); % Creating connectivity matrix for CA3-CA3 synapses
AC_Connectivity_Ratio = AC_Connectivity_Ratio - diag(diag(AC_Connectivity_Ratio)); % Setting the connectivity probabillity of a neuron to itself to zero

neuronGroup.AC_Connectivity_Ratio = AC_Connectivity_Ratio;
% tpulse=0.1;       % Time when pulse starts
% lengthpulse=0.05;
Nt=settings.tmax/settings.dt;       % Total number of time-steps. 



% [spikeMat, tVec] = ISI_IBI_spikeGen(tmax, num_of_presynaptic_neurons, dt , muA, muB, sigmaA, sigmaB, fr);

% This is the simplest simulation with a poisson spike generator. Each
% neuron will fire with a poisson firing rate defined by the firing
% frequency.
if ~settings.one_on_one 
    [spikeMat, tVec,h] = poissonSpikeGen(params, settings);
end


% A simulation where there is one presynaptic neuron and one post synaptic
% neuron. The spikeMat generated will be 1xNt.
% The presynaptic neuron will fire 7 times total.
% The first 6 times will be with a time interval of 200 time steps between them.
% The 7th time will be much later. This way we get the behaviour that a
% very long time has passed since the intense simulation.

if settings.one_on_one;
    spikeMat = zeros(1,Nt);
    DT = 200; %
    for i = DT : DT : 6*DT 
        spikeMat(i) = 1;
    end
    tVec = (0:settings.dt:settings.tmax-settings.dt)*1e3;
    spikeMat(40*DT) = 1;
    h = plotRaster(logical(spikeMat),tVec);
    
end

% This is generalization of the one_on_one simulation. Here we have several
% presynaptic neurons firing with the same firing pattern as was defined in
% the one_on_one state. 
if settings.one_on_one_multi
    spikeMat = zeros(1,Nt);
    DT = 200; %
    for i = DT : DT : 6*DT 
        spikeMat(i) = 1;
    end
    tVec = (0:settings.dt:settings.tmax-settings.dt)*1e3;
    spikeMat(40*DT) = 1;
    spikeMat = repmat(spikeMat,params.num_of_presynaptic_neurons,1);
    h = plotRaster(logical(spikeMat),tVec);
end

% A simulation with account for Inter-Spike-Intervals and
% Inter-Burst-Intervals.
if settings.ISI_IBI
    [spikeMat, tVec,h] = ISI_IBI_spikeGen(params, settings);
end
    
    
saveas(h,'Figures\Pre Synaptic Raster Plot.bmp');



neuronGroup.Pre_Synaptic_Spikes = spikeMat;
neuronGroup.Pre_Synaptic_Spikes_Times = tVec;
% for i = 1 : round(lengthpulse/dt)
%     temp_mat = [zeros(size(spikeMat,1),1) spikeMat(:,1:end-1)];
%     spikeMat = spikeMat + temp_mat;
% end
% tempI = logical(spikeMat) * Iapp;




% tempI is the synaptic current casued by each presynaptic neuron
% individually. I will be the current felt by each postsynaptic neuron
% which is the linear sum of all presynaptic neurons connected to it.





% We used two version to solve the differential equations regarding u and x.
% The first form was derived from the paper: Differential signaling via the same axon of neocortical
% pyramidal neurons. In this form, there is no parameter f.
% The second form was taken from the paper Probabilistic inference of short-term synaptic plasticity in neocortical microcircuits, Costa et al., 2013
% Please be aware of this for you will have to change which function of
% solving is used.

% Building the EPSC matrix using the model
% Version 1
% % tempI = solve_u_x (spikeMat,tmax,dt,A,U,tauD,tauF);

% Version 2
[tempI,u,x] = solve_u_x_Meish(spikeMat,tVec, params, settings);
neuronGroup.u = u;
neuronGroup.x = x;



for i = 1 : size(Connectivity_Ratio,2)
    idx = find(Connectivity_Ratio(:,i));
    curI = sum(tempI(idx,:),1);
    try
        I = [I; curI];
    catch
        I = curI;
    end
        
end
I(I<0) = 0; % Making sure there is no IPSC




T=0:settings.dt:settings.tmax; 
% I=zeros(num_of_neurons,length(T));
gsra=zeros(params.num_of_postsynaptic_neurons,length(T));    % spike-rate adaptaion conductance
tau_sra = 0.2;          % time for spike-rate adaptation to decay
delta_gsra = 5e-9;     % increase in spike-rate conductance per spike

gref=zeros(size(T));    % spike refractory conductance
tau_ref = 0.002;        % time for refractory conductance to decay
delta_gref = 200e-9;     % increase in refractory conductance per spike

V = ones(params.num_of_postsynaptic_neurons,1) * params.E;                 % begin simulation with membrane potential at leak level
gsra = zeros(params.num_of_postsynaptic_neurons,1);            % initially no spike-rate adaptation
gref = zeros(params.num_of_postsynaptic_neurons,1);            % initially no refractory current

spikes=zeros(params.num_of_postsynaptic_neurons,Nt);  % vector to record spike times

% fprintf('A current with amplitude of %i nA is injected from %i to %i \n',Iapp*1e9, (tpulse/dt), (tpulse+lengthpulse)/dt)


% I(:,round(tpulse/dt) : round((tpulse+lengthpulse)/dt)) = Iapp;


IAC = zeros(params.num_of_postsynaptic_neurons,1);
neuronGroup.Post_Synaptic_Spike_Times = zeros(params.num_of_postsynaptic_neurons,Nt);
for i = 2:Nt
    
    
% %     Calculating the synaptic current felt by the AC connections

    
    IAC(:,end) = IAC(:,end) * params.IACsize;
    
    % Next line solves tau_sra*(d gsra/dt) = -gsra
    gsra(:,i) = gsra(:,i-1)*(1-settings.dt/tau_sra);
    % Next line solves tau_ref*(d gref/dt) = -gref
    gref(:,i) = gref(:,i-1)*(1-settings.dt/tau_ref);
    
    
    V = [V V(:,i-1) + settings.dt/params.tau_m*(params.E-V(:,i-1) + params.Rm * I(:,i-1) -...
        (V(:,i-1)-params.E_K)*params.Rm.*(gsra(:,i-1)+gref(:,i-1))) + params.Rm * IAC(:,i-1)];
    V(V(:,i-1) == params.Vspike,i) = params.Vreset;
    
    
    thresh_idx = V(:,i) > params.Vth; % Finding the neurons for which the potential exceeds the threshold potential
    neuronGroup.Post_Synaptic_Spike_Times(thresh_idx,i) = 1;
    
    V(thresh_idx,i) = params.Vspike; 
%     spikes(thresh_idx,i) = 1; 
    tmpIAC = zeros(params.num_of_postsynaptic_neurons,1);
    tmpIAC(thresh_idx) = 1; % An AC current is genrated at the exact time of the neurons AP
    try 
        IAC = [IAC AC_Connectivity_Ratio * tmpIAC];
    catch
        IAC = AC_Connectivity_Ratio * tmpIAC;
    end
%     The refractory period is created as a period in which the current
%     reaching the neuron is set to zero. This way we still get the
%     behaviour of the neuron following the AP but we will have a
%     refractory period.

    I(thresh_idx,i:i+round(params.refract*1e-3/settings.dt)) = 0;
    IAC(thresh_idx,i:round(i+params.refract*1e-3/settings.dt)) = 0;
    
    gsra(thresh_idx,i) = gsra(thresh_idx,i-1) + delta_gsra; 
    gref(thresh_idx,i) = gref(thresh_idx,i-1) + delta_gref;
    
    
end
spikes(I~=0) = 1;
neuronGroup.Spike_Times = spikes;
neuronGroup.Post_Synaptic_Potential = V;
neuronGroup.Post_Synaptic_Current = I;


