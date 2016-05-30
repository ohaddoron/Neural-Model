close all; clear all; clc;
tSim = 1;
nTrials = 20;
dt = 0.0002;
muA = 25e-3;
muB = 150e-3;
sigmaA = 2e-2;
sigmaB = 1.5*sigmaA;
fr = 25;
num_of_tests = 1;
m = nan(num_of_tests,1);
A = 2;

binWidth = 50*dt;
bin = 0:binWidth:(tSim-dt);
val = zeros(1,length(bin));

for i = 1 : num_of_tests
    
    
    [spikeMat, tVec] = ISI_IBI_spikeGen(tSim, nTrials, dt , muA, muB, sigmaA, sigmaB, fr, A);
%     close all;
    m(i) = round(mean(sum(spikeMat,2)),1);
    disp(num2str(i))
    disp(num2str(m(i)));
    for j = 1 : size(spikeMat,1)
        d = diff(find(spikeMat(j,:)));
        Y = discretize(d*dt,bin);
        val(Y) = val(Y) + 1;
    end
%     close all;
end
figure;
hist(m);



figure;

bar(bin*1e3,val);
axis tight;

    