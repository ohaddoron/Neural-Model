stepSize = settings.SW/settings.dt;
b = 0;
for n = 1 : (numel(tVec) - stepSize)
    spikes = a.Post_Synaptic_Spike_Times(:,n:n+stepSize-1);
    b = b + sum(spikes(:));
end