function [EPSC,u_mat,x_mat] = solve_u_x_Meish (spikeMat,tVec, params, settings)

% In this version the equations used to solve script were taken from the 
% paper Probabilistic inference of short-term synaptic plasticity in neocortical microcircuits, Costa et al., 2013
% In this version another parameter is added which is f


Nt = settings.tmax/settings.dt; % the number of iterations in the time domain.

% The solution is applied to each neuron individually. The problem is that
% the AP for each neuron occur at different times so the DeltaT for each of
% them is different. It's a problem which will be worth solving in the
% future for simulations that will need to run for a longer time.

x_mat = zeros(size(spikeMat,1), size(spikeMat,2));
u_mat = zeros(size(spikeMat,1), size(spikeMat,2));
for i = 1 : size(spikeMat,1)
    
    spikeTrain = spikeMat(i,:); % Getting the AP for the current neuron
    APtimes = find(spikeTrain); % finding the last time that an AP occured
    if (length(APtimes) >=2) 
        DeltaT = (APtimes(2:end) - APtimes(1:end-1))*settings.dt;
    elseif (length(APtimes) == 1) 
%         fprintf('only 1 AP in pre-synaptic neuron number %i', i);
    else 
%         fprintf('pre-synaptic neuron number %i did not fire APs', i);
        continue; 
    end
    u_vec = zeros(1, length(APtimes));
    x_vec = zeros(1, length(APtimes));
    % actually we can even vectorize this part to make it even more
    % efficient
        for ap= 1 : length(APtimes)
            if ap ~= 1;
                 u_vec(ap) = params.U + (u_vec(ap-1) + params.f*(1-u_vec(ap-1)) - params.U)*exp(-DeltaT(ap-1)/params.tauF);
                 x_vec(ap) = 1- (1 - x_vec(ap-1) * (1- u_vec(ap-1))) * exp(-DeltaT(ap-1)/params.tauD);
            else
                u_vec(ap) = params.U;
                x_vec(ap) = 1;
            end
        end
        assert(length(u_vec) == length(APtimes));
        % insert the x and u values to their corresponding matrices
        % according to their time of occurence.
        x_mat(i,APtimes) = x_vec; 
        u_mat(i,APtimes) = u_vec;
    
end


h = figure('Visible','Off');
if size(x_mat,1) > 4 % if the number of pre-synaptic neurons is more than 4
    for i = 1 : 4
        subplot(2,2,i)
        hold on;
        % plot for each of them
        tmp_x_mat = x_mat(i,:);
        tmp_u_mat = u_mat(i,:);
        tmp_tVec = tVec;
        tmp_tVec(x_mat(i,:) == 0) = [];
        tmp_x_mat(x_mat(i,:) == 0) = [];
        tmp_u_mat(x_mat(i,:) == 0) = [];
        plot(tmp_tVec,tmp_x_mat,'-o');  % x
        plot(tmp_tVec,tmp_u_mat,'-o');  % u 
        plot(tmp_tVec,params.Rm*params.A*tmp_x_mat.*tmp_u_mat,'-o'); % and x*u*A
    %     stem(spikeMat(i,:),'.','k'); % This will mark the AP
        if i == 1
            legend('R','u','EPSC','Location','best');
        end
        set(gca,'XLim',[0 settings.tmax*1e3]);

    end
else
    hold on;
    plot(tVec,x_mat(1,:));
    plot(tVec,u_mat(1,:));
    plot(tVec,x_mat(1,:).*u_mat(1,:).*params.A);
    legend('R','u','EPSC','Location','best');
    axis tight;
end

saveas(h,'Figures\u_x.bmp');
EPSC = [];
for i = 1 : size(spikeMat,1)
    tmp = awgn(params.A*x_mat(i,:).*u_mat(i,:),10,'measured');
    if ~settings.one_on_one && ~settings.one_on_one_multi
        EPSC = [EPSC; spikeMat(i,:).*(tmp)];
    else
        EPSC = [EPSC; spikeMat(i,:).*(params.A*(x_mat(i,:).*u_mat(i,:)))];
    end
end


