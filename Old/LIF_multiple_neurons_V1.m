% Leaky-integrate-and-fire neuron with spike rate adaptation and a
% refractory current.

clear all; close all; clc;

tic

num_of_neurons = 500;

E=-70e-3;
Rm=10e6;
tau_m=10e-3;
Vth=-54e-3; 
Vreset=-80e-3;
Vspike = 20e-3;
E_K = -80e-3;       % reversal potential for spike-rate adaptation (SRA) and 
                    % refractory current

Iapp = 4e-9;      % amplitude of injected current  
dt=0.0002;        
tmax=0.7;
tpulse=0.1;       % Time whe pulse starts
lengthpulse=0.5;
Nt=tmax/dt;       % Total number of time-steps. 
T=0:dt:tmax; 
I=zeros(num_of_neurons,length(T));
gsra=zeros(num_of_neurons,length(T));    % spike-rate adaptaion conductance
tau_sra = 0.2;          % time for spike-rate adaptation to decay
delta_gsra = 5e-9;     % increase in spike-rate conductance per spike

gref=zeros(size(T));    % spike refractory conductance
tau_ref = 0.002;        % time for refractory conductance to decay
delta_gref = 200e-9;     % increase in refractory conductance per spike

V = ones(num_of_neurons,1) * E;                 % begin simulation with membrane potential at leak level
gsra = zeros(num_of_neurons,1);            % initially no spike-rate adaptation
gref = zeros(num_of_neurons,1);            % initially no refractory current

spikes=zeros(num_of_neurons,length(T));  % vector to record spike times

fprintf('A current with amplitude of %i nA is injected from %i to %i \n',Iapp*1e9, (tpulse/dt), (tpulse+lengthpulse)/dt)


I(:,round(tpulse/dt) : round((tpulse+lengthpulse)/dt)) = Iapp;

disp('Nt:');
disp(Nt);

for i = 2:Nt;
    
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
toc
%clf;
% figure(1);
% subplot(3,1,1);
% plot(T,I); % Current as a function of time.
% ylabel('Current [A]');
% hold on;
% subplot(3,1,2);  %Voltage as a function of time.
% plot(T,V);
% ylabel('Voltage');
% hold on;
% subplot(3,1,3);
% plot(T,gsra);    %spike-rate adpatation conductance as a function of time
% hold on
% plot(T,gref,'r'); %Refractory conductance as a function of time
% legend('gsra', 'gref');