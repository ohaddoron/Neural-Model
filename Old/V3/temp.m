clear all; close all; clc;
tic

tmax = 3;
dt = 0.0002;
Nt = tmax/dt;

[spikeMat, tVec] = poissonSpikeGen(10, tmax, 500,dt);
tmp = double(spikeMat);
tmp(tmp == 0) = nan;
lastAP = nan(size(spikeMat,1),1);
% DeltaT = nan(size(spikeMat));
DeltaT = [];
for n = 1 : Nt
    [M,idx] = max(fliplr(tmp(:,1:n)),[],2);
    idx(isnan(M)) = nan;
    %nansum(~isnan(idx))
    idx = size(tmp(:,1:n),2) - idx + 1;
    curDeltaT = idx - lastAP;
    curDeltaT(isnan(curDeltaT)) = idx(isnan(curDeltaT));
    curDeltaT = curDeltaT * dt;
    DeltaT = [DeltaT curDeltaT];
    lastAP(~isnan(curDeltaT)) = curDeltaT(~isnan(curDeltaT));
    
    
    
end
    
toc