function [spikeMat, tVec] = ISI_IBI_spikeGen(tSim, nTrials, dt , muA, muB, sigmaA, sigmaB,p , fr)

% This function will create a presynaptic spikes matrix which is dependent
% on inter-spike-interval and inter-burst-interval times. 
% This function will draw a number from a bernoulli distribution with the probability
% of success being dependent on the time passed since the last spike.
% tSim is the simulation time in seconds.
% dt is the size of each time step in the simulation.
% nTrials is the number of presynaptic neurons.
% p is the default chance for a spike


nBins = floor(tSim/dt);
spikeMat = zeros(nTrials,nBins);
tVec = (0:dt:tSim-dt)/dt;


X = 0:dt:tSim;
Y = normpdf(X,muA,sigmaA) + normpdf(X,muB,sigmaB); 
% Y is the sum of two normal distributions with different means and std.
% The values in Y are probabilities of spikes as a function of time
% intervals.
figure;
hold on;
Y = Y/(sum(Y));
plot(Y);
B = cumprod(1-Y)./(1-Y) .* Y;
plot(B);
Y = 4 * fr * Y * sum(X.*B);
plot(Y);
figure;
area(X,Y);


for i = 1 : size(spikeMat,1)
    for j = 1 : size(spikeMat,2)
        lastAP = find(spikeMat(i,1:j-1),1,'last');
        if ~isempty(lastAP)
            P = Y(j-lastAP); % The probability P is the value in Y which corresponds to how long has it been since the last spike
            spikeMat(i,j) = binornd(1,P);
        else
            spikeMat(i,j) = binornd(1,p);
        end
    end
end
spikeMat = logical(spikeMat);
figure;
imagesc(~spikeMat);
colormap(gray);

end

