function [spikeMat, tVec,h] = ISI_IBI_spikeGen(params,settings)

% This function will create a presynaptic spikes matrix which is dependent
% on inter-spike-interval and inter-burst-interval times. 
% This function will draw a number from a bernoulli distribution with the probability
% of success being dependent on the time passed since the last spike.
% tSim is the simulation time in seconds.
% dt is the size of each time step in the simulation.
% nTrials is the number of presynaptic neurons.
% p is the default chance for a spike
A = 1.85;

nBins = floor(settings.tmax/settings.dt);
spikeMat = zeros(params.num_of_presynaptic_neurons,nBins);
tVec = (0:settings.dt:settings.tmax-settings.dt);


X = 0:settings.dt:settings.tmax;
Y = normpdf(X,params.muISI,params.sigmaISI) + normpdf(X,params.muIBI,params.sigmaIBI); 
% Y is the sum of two normal distributions with different means and std.
% The values in Y are probabilities of spikes as a function of time
% intervals.
Y = Y/(sum(Y));
B = cumprod(1-Y)./(1-Y) .* Y;
Y = 1/(sum(X.*B)) * B;

Y = A * Y/max(Y);
% figure;
% area(X*1e3,Y*settings.dt*params.fr);

for i = 1 : size(spikeMat,1)
    for j = 1 : size(spikeMat,2)
        lastAP = find(spikeMat(i,1:j-1),1,'last');
% if the neuron has already fired before, the next time it fires will be by
% the the time passed since the last fire, otherwise, the lastAP variable
% will be drawn from a normal distribution with mean of 200 and std of 100.
        if ~isempty(lastAP) 
            P = Y(j-lastAP); % The probability P is the value in Y which corresponds to how long has it been since the last spike
        else
            
            if j == 1
                pseudo_lastAP = inf;
                while (j-pseudo_lastAP) < 1 || (j-pseudo_lastAP) > numel(Y)
                    pseudo_lastAP =  -round(normrnd(params.mu_lastAP/settings.dt * 1e-3,params.sigma_lastAP/settings.dt * 1e-3)); % last time that the neuron fired drawn from a normal distribution
                end
            end
            if (j-pseudo_lastAP) < numel(Y)
                P = Y(j - pseudo_lastAP);
            else
                P = 0;
            end
        end
            spikeMat(i,j) = rand(1,1) <     P * params.fr * settings.dt;
%             spikeMat(i,j) = rand(1, 1) < params.fr*settings.dt; % This
%             was used before the pseudo_lastAP parameters was entered
    end
end
spikeMat = logical(spikeMat);
h=plotRaster(spikeMat,tVec*1e3);
% h.Visible = 'on';
% saveas(h,'Figures\Raster Plot.bmp');
% figure;
% imagesc(~spikeMat);
% set(gca,'XTickLabel',round(dt*tVec(get(gca,'XTick')),2));
% colormap(gray);

end

