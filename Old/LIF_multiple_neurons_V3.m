% Leaky-integrate-and-fire neuron with spike rate adaptation and a
% refractory current.

clear all; close all; clc;
tic

num_of_presynaptic_neurons = 500;
num_of_postsynaptic_neurons = 100;
p = 0.01;

Connectivity_Ratio = connect(num_of_presynaptic_neurons, num_of_postsynaptic_neurons, p);



E=-70e-3;
Rm=10e6;
tau_m=10e-3;
Vth=-54e-3; 
Vreset=-80e-3;
Vspike = 20e-3;
E_K = -80e-3;       % reversal potential for spike-rate adaptation (SRA) and 
                    % refractory current

Iapp = 4e-10;      % size of each synaptic current. This is an arbitrary size.
dt=0.0002;        
tmax=1;
% tpulse=0.1;       % Time when pulse starts
lengthpulse=0.05;
Nt=tmax/dt;       % Total number of time-steps. 

[spikeMat, tVec] = poissonSpikeGen(10, tmax, num_of_presynaptic_neurons,dt);


% for i = 1 : round(lengthpulse/dt)
%     temp_mat = [zeros(size(spikeMat,1),1) spikeMat(:,1:end-1)];
%     spikeMat = spikeMat + temp_mat;
% end
% tempI = logical(spikeMat) * Iapp;

% Building the EPSC matrix using the model


% tempI is the synaptic current casued by each presynaptic neuron
% individually. I will be the current felt by each postsynaptic neuron
% which is the linear sum of all presynaptic neurons connected to it.



    
tempI = solve_u_x (spikeMat,tmax,dt);

for i = 1 : size(Connectivity_Ratio,2)
    idx = find(Connectivity_Ratio(:,i));
    curI = sum(tempI(idx,:),1);
    try
        I = [I; curI];
    catch
        I = curI;
    end
        
end

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

spikes=zeros(num_of_postsynaptic_neurons,length(T));  % vector to record spike times

% fprintf('A current with amplitude of %i nA is injected from %i to %i \n',Iapp*1e9, (tpulse/dt), (tpulse+lengthpulse)/dt)


% I(:,round(tpulse/dt) : round((tpulse+lengthpulse)/dt)) = Iapp;

disp('Nt:');
disp(Nt);

for i = 2:Nt
    
    % Next line solves tau_sra*(d gsra/dt) = -gsra
    gsra(:,i) = gsra(:,i-1)*(1-dt/tau_sra);
    % Next line solves tau_ref*(d gref/dt) = -gref
    gref(:,i) = gref(:,i-1)*(1-dt/tau_ref);
    
    
    V = [V V(:,i-1) + dt/tau_m*(E-V(:,i-1) + Rm * I(:,i-1) -...
        (V(:,i-1)-E_K)*Rm.*(gsra(:,i-1)+gref(:,i-1)))];
    V(V(:,i-1) == Vspike,i) = Vreset;
    
    
    thresh_idx = V(:,i) > Vth; % Finding the neurons for which the potential exceeds the threshold potential
    
    V(thresh_idx,i) = Vspike; 
    spikes(thresh_idx,i) = 1; 
    gsra(thresh_idx,i) = gsra(thresh_idx,i-1) + delta_gsra; 
    gref(thresh_idx,i) = gref(thresh_idx,i-1) + delta_gref;
    
    
end
% Plotting the post synaptic neurons potential and synaptic current
figure;
suptitle('Post Synaptic Potential')
for i = 1 : 16
    subplot(4,4,i)
    plot(V(i,:)*1e3);
    axis tight
    ylabel('Membrane potential [mV]');
    ax = gca;
    ax.XTickLabel = num2cell(str2num(cell2mat(ax.XTickLabel)) * dt );
    xlabel('Time [sec]');
end
figure
suptitle('Synaptic Current')
for i = 1 : 16
    subplot(4,4,i)
    plot(I(i,:)*1e9);
    ylabel('Synaptic current [nA]');
    set(gca,'YLim',[0 max(I(:)*1e9)]);
    set(gca,'XTickLabel',ax.XTickLabel);
    xlabel('Time [sec]');
end

toc