function [h] = plotRaster(spikeMat, tVec)

% This functuon will create a raster plot for the neurons. 
% spikeMat is a logical matrix with 1 at the time point which a spike
% occurs and 0 otherwise. 
% tVec is the time vector corresponding to each column in spikeMat
% The code was taken from https://praneethnamburi.wordpress.com/2015/02/05/simulating-neural-spike-trains/
h = figure('Visible','Off');
tVec = tVec * 1e3;
hold all;
for trialCount = 1:size(spikeMat,1)
    spikePos = tVec(spikeMat(trialCount, :));
    for spikeCount = 1:length(spikePos)
        plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.4 trialCount+0.4], 'k');
    end
end
ylim([0 size(spikeMat, 1)+1]);
title('Raster plot');
set(gca,'XTickLabel',tVec);
xlabel('Time (ms)');
ylabel('Trial Number');