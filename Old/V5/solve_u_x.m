function [EPSC,u,x] = solve_u_x (spikeMat,tmax,dt,A,U,tauD,tauF,f,tVec)

% In this version the equations used to solve script were taken from the 
% paper Probabilistic inference of short-term synaptic plasticity in neocortical microcircuits, Costa et al., 2013
% In this version another parameter is added which is f


Nt = tmax/dt;

% Creating the AP matrix



u = [];
x = [];
% The solution is applied to each neuron individually. The problem is that
% the AP for each neuron occur at different times so the DeltaT for each of
% them is different. It's a problem which will be worth solving in the
% future for simulations that will need to run for a longer time.

for i = 1 : size(spikeMat,1)
    
    spikeTrain = spikeMat(i,:); % Getting the AP for the current neuron
    cur_u(1) = U; % Initializing the synaptic parameters
    cur_x(1) = 1-U;

    for n = 1 : Nt-1
        % Calculating the time difference between the nth AP and the
        % (n+1)th AP.
%         APtimes = find(spikeTrain(1:n),2,'last');
%         if length(APtimes) < 1
%             DeltaT = inf;
%         else
%             if length(APtimes) == 2
%                 DeltaT = (APtimes(2)-APtimes(1))*dt;
%             else
%                 DeltaT = APtimes * dt;
%             end
%         end
        APtimes = find(spikeTrain(1:n),1,'last'); % finding the last time that an AP occured
        if ~isempty(APtimes)
            DeltaT = (n - APtimes)*dt; % This is the time difference (in seconds) between the current time and the last AP
        else
            DeltaT = inf; % if there hasn't been an AP yet, the time difference is infinity.
        end
        % Solving the equations for u and x.
        cur_u(n+1) = U + (cur_u(n) + f*(1-cur_u(n)) - U)*exp(-DeltaT/tauF);
        cur_x(n+1) = 1- (1 - cur_x(n) * (1- cur_u(n))) * exp(-DeltaT/tauD);
    end
    u = [u; cur_u];
    x = [x; cur_x];
    
end


% The model is currently working for 9 neurons so we can visualize it
% properly. 
h = figure('Visible','Off');
if size(x,1) > 4
    for i = 1 : 4
        subplot(2,2,i)
        hold on;

        plot(tVec,x(i,:));
        plot(tVec,u(i,:));
        plot(tVec,A*x(i,:).*u(i,:));
    %     stem(spikeMat(i,:),'.','k'); % This will mark the AP
        if i == 1
            legend('R','u','EPSC','Location','northwestoutside');
        end
        axis tight

    end
else
    hold on;
    plot(tVec,x(1,:));
    plot(tVec,u(1,:));
    plot(tVec,x(1,:).*u(1,:));
    legend('R','u','EPSC','Location','best');
    axis tight;
end

saveas(h,'Figures\u_x.bmp');
EPSC = [];
for i = 1 : size(spikeMat,1)
    tmp = awgn(A*x(i,:).*u(i,:),10,'measured');
%     EPSC = [EPSC; spikeMat(i,:).*(tmp)];
    EPSC = [EPSC; spikeMat(i,:).*(A*(x(i,:).*u(i,:)))];
end









    


    