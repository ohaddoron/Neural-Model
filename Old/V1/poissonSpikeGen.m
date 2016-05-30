function [spikeMat, tVec] = poissonSpikeGen(fr, tSim, nTrials, dt)

% This function will create the presynatpic firing rate. 
% The presynaptic neurons fire with a poisson disturbiuton.
% fr - The firing rate
% tSim - Length of simulation in seconds
% nTrials - number of presynaptic neurons
% dt - the time interval used
% % Notice that dt was defined here as 1ms, we may want to change this to fit
% % our simulation
% spikeMat is a 2D matrix, with the size of nTrails x tSim/dt. Each column
% corresponds to a time point and each row is a presynatpic neuron.
% tVec is a time vector describing the time stamp for each column in
% spikeMat.
% The code was taken from https://praneethnamburi.wordpress.com/2015/02/05/simulating-neural-spike-trains/


nBins = floor(tSim/dt);
spikeMat = rand(nTrials, nBins) < fr*dt;
tVec = (0:dt:tSim-dt)/dt;
h = figure;
imagesc(~spikeMat);
colormap(gray);
title('Raster Plot');
saveas(h,'Figures\Raster Plot.bmp');

%plotRaster(spikeMat,tVec);