function neuronGroup = simulate(num_of_presynaptic_neurons,num_of_postsynaptic_neurons,...
    p,fr,E,Rm,tau_m,Vth,Vreset,Vspike,E_K,dt,tmax,A,U,tauD,tauF, pAC, IACsize, f)

% Leaky-integrate-and-fire neuron with spike rate adaptation and a
% refractory current.





Connectivity_Ratio = connect(num_of_presynaptic_neurons, num_of_postsynaptic_neurons, p);
neuronGroup.Connectivity_Ratio = Connectivity_Ratio;
neuronGroup.num_of_presynaptic_neurons = num_of_presynaptic_neurons;
neuronGroup.num_of_postsynaptic_neurons = num_of_postsynaptic_neurons;

AC_Connectivity_Ratio = connect(num_of_postsynaptic_neurons, num_of_postsynaptic_neurons, pAC); % Creating connectivity matrix for CA3-CA3 synapses
AC_Connectivity_Ratio = AC_Connectivity_Ratio - diag(diag(AC_Connectivity_Ratio)); % Setting the connectivity probabillity of a neuron to itself to zero

neuronGroup.AC_Connectivity_Ratio = AC_Connectivity_Ratio;
% tpulse=0.1;       % Time when pulse starts
% lengthpulse=0.05;
Nt=tmax/dt;       % Total number of time-steps. 


muISI = 25e-3;
muB = 250e-3;
sigmaA = 3e-2;
sigmaB = 1.5 * sigmaA;

% [spikeMat, tVec] = ISI_IBI_spikeGen(tmax, num_of_presynaptic_neurons, dt , muA, muB, sigmaA, sigmaB, fr);
[spikeMat, tVec,h] = poissonSpikeGen(fr, tmax, num_of_presynaptic_neurons, dt);

spikeMat = zeros(1,length(spikeMat));
DT = 100;
for i = 1 : DT : 6*DT 
    spikeMat(i) = 1;
end
tVec = (0:dt:tmax-dt)*1e3;
spikeMat(20*DT) = 1;
h = plotRaster(logical(spikeMat),tVec);
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
[tempI,u,x] = solve_u_x (spikeMat,tmax,dt,A,U,tauD,tauF , f,tVec);
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




T=0:dt:tmax; 
% I=zeros(num_of_neurons,length(T));
gsra=zeros(num_of_postsynaptic_neurons,length(T));    % spike-rate adaptaion conductance
tau_sra = 0.2;          % time for spike-rate adaptation to decay
delta_gsra = 5e-9;     % increase in spike-rate conductance per spike

gref=zeros(size(T));    % spike refractory conductance
tau_ref = 0.002;        % time for refractory conductance to decay
delta_gref = 200e-9;     % increase in refractory conductance per spike

V = ones(num_of_postsynaptic_neurons,1) * E;                 % begin simulation with membrane potential at leak level
gsra = zeros(num_of_postsynaptic_neurons,1);            % initially no spike-rate adaptation
gref = zeros(num_of_postsynaptic_neurons,1);            % initially no refractory current

spikes=zeros(num_of_postsynaptic_neurons,Nt);  % vector to record spike times

% fprintf('A current with amplitude of %i nA is injected from %i to %i \n',Iapp*1e9, (tpulse/dt), (tpulse+lengthpulse)/dt)


% I(:,round(tpulse/dt) : round((tpulse+lengthpulse)/dt)) = Iapp;

disp('Nt:');
disp(Nt);
IAC = zeros(num_of_postsynaptic_neurons,1);
for i = 2:Nt
    
    
% %     Calculating the synaptic current felt by the AC connections

    
    IAC(:,end) = IAC(:,end) * IACsize;
    
    % Next line solves tau_sra*(d gsra/dt) = -gsra
    gsra(:,i) = gsra(:,i-1)*(1-dt/tau_sra);
    % Next line solves tau_ref*(d gref/dt) = -gref
    gref(:,i) = gref(:,i-1)*(1-dt/tau_ref);
    
    
    V = [V V(:,i-1) + dt/tau_m*(E-V(:,i-1) + Rm * I(:,i-1) -...
        (V(:,i-1)-E_K)*Rm.*(gsra(:,i-1)+gref(:,i-1))) + IAC(:,i-1)];
    V(V(:,i-1) == Vspike,i) = Vreset;
    
    
    thresh_idx = V(:,i) > Vth; % Finding the neurons for which the potential exceeds the threshold potential
    
    V(thresh_idx,i) = Vspike; 
%     spikes(thresh_idx,i) = 1; 
    tmpIAC = zeros(num_of_postsynaptic_neurons,1);
    tmpIAC(thresh_idx) = 1; % An AC current is genrated at the exact time of the neurons AP
    try 
        IAC = [IAC AC_Connectivity_Ratio * tmpIAC];
    catch
        IAC = AC_Connectivity_Ratio * tmpIAC;
    end
    gsra(thresh_idx,i) = gsra(thresh_idx,i-1) + delta_gsra; 
    gref(thresh_idx,i) = gref(thresh_idx,i-1) + delta_gref;
    
    
end
spikes(I~=0) = 1;
neuronGroup.Spike_Times = spikes;
neuronGroup.Post_Synaptic_Potential = V;
neuronGroup.Post_Synaptic_Current = I;
% Plotting the post synaptic neurons potential and synaptic current


