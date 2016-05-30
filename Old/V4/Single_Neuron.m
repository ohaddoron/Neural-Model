function neuron_group = Single_Neuron (spikes,num_of_neurons,run_time) 

% This function will solved the differential equations which describe the
% neuron's behaviour over time.
% spikes is a vector with the EPSP caused by a presynaptic spike
% The neural model is taken from Adaptive Exponential Integrate-and-Fire Model as an Effective Description of Neuronal Activity, J. Neurophysiol., 2005

%%
C = 281; % pF
gL = 30; % nS
taum = C/gL; % msec
EL = -70.6; % mV
VT = -50.4; % mV
DeltaT = 2; % mV
Vcut = VT + 5*DeltaT; % mV
tauw = 144; % msec
a = 4; % nS
b = 0.0805; % nA
Vr  = -70.6; % mV
dt = 1e-3;
I = 0;

syms V(t) w(t)

eqn1 = diff(V,t) == (gL*(EL-V) + gL * DeltaT * exp((V-VT)/DeltaT) + I - w)/C;
eqn2 = diff(w,t) == (a*(V-EL)-w)/tauw;

S = dsolve(eqn1,eqn2);
VSol(t) = S.V
wSol(t) = S.w
