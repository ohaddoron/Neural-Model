function [Sync_Indicator,num_of_post_synaptic_spikes,sub_treshold_EPSP] = plotFigures (neuronGroup,params,settings)
    Sync_Indicator = zeros(numel(neuronGroup),1);
    num_of_post_synaptic_spikes = zeros(1,numel(neuronGroup));
    sub_treshold_EPSP = zeros(1,numel(neuronGroup));
    tVec = 0:settings.dt:settings.tmax-settings.dt;
    
    for i = 1 : numel(neuronGroup)
        if numel(params.pAC) > 1
            outPath = fullfile(pwd,'Figures',strcat('pAC',num2str(params.pAC(i))));
        else
            switch i
                case 1
                    path2 = fullfile(pwd,'Strong Depression');
                case 2
                    path2 = fullfile(pwd,'Depression');
                case 3
                    path2 = fullfile(pwd,'Facilitation-Depression');
                case 4
                    path2 = fullfile(pwd,'Facilitation');
                case 5
                    path2 = fullfile(pwd,'Strong Facilitation');
            end
            outPath = path2;
        end
        mkdir(outPath);
        
        num_of_post_synaptic_spikes(i) = sum(neuronGroup{i}.Post_Synaptic_Spike_Times(:));
        
            
        
        V = neuronGroup{i}.Post_Synaptic_Potential;
        I = neuronGroup{i}.Post_Synaptic_Current;
        h = figure('Visible','Off');
        suptitle('Post Synaptic Potential')
        idx = find(mean(abs(V-params.E),2) > 0, 4);
        if numel(idx) < 4
            idx =1:4;
        end
        if size(V,1) >= 4
            num_of_plots = 4;
        else
            num_of_plots = length(idx);
        end
        for ii = 1 : num_of_plots
            if num_of_plots < 4
                subplot(1,1,ii)
            else
                subplot(2,2,ii)
            end
            plot((0:length(V(idx(ii),:))-1)*settings.dt*1e3,V(idx(ii),:)*1e3);
            axis tight
            ylabel('Potential[mV]');
            ax = gca;
        %     ax.XTickLabel = num2cell(str2num(cell2mat(ax.XTickLabel)) * dt );
            xlabel('Time [msec]');
        end

        

        saveas(h,fullfile(outPath,'Post Synaptic Potential.bmp'));

        h = figure('Visible','Off');
        suptitle('Synaptic Current')
        if size(I,1) > 4
            for ii = 1 : 4
                subplot(2,2,ii)
                plot((0:length(I(idx(ii),:))-1)*settings.dt*1e3,I(idx(ii),:)*1e9);
                ylabel('Current [nA]');
                set(gca,'YLim',[0 max(I(:)*1e9)]);
            %     set(gca,'XTickLabel',ax.XTickLabel);
                xlabel('Time [msec]');
            end
        else
            
                plot((0:length(I(1,:))-1)*settings.dt*1e3,I(1,:)*1e9);
                ylabel('Current [nA]');
                set(gca,'YLim',[0 max(I(:)*1e9)]);
            %     set(gca,'XTickLabel',ax.XTickLabel);
            xlabel('Time [msec]');
        end
        saveas(h,fullfile(outPath,'Synaptic Current.bmp'));
        
        
        V = neuronGroup{i}.Post_Synaptic_Potential - params.E;
        
        
        
        
        
        V(V == params.Vspike-params.E) = nan;

        h = figure('Visible','Off');
        hold all;
        binWidth = 50;
        bin = 0:binWidth:settings.tmax*1e3;
        val = nan(1,length(bin));
        h = figure('Visible','Off');
        hold all;
        for m = 1 : size(V,1)
            [pks,locs] = findpeaks(V(m,:)*1e3); % This will locate the location and size of the EPSPs.
            if ~isempty(pks)
                sub_treshold_EPSP(i) = sum([sub_treshold_EPSP(i),sum(pks(:))]);
                locs = locs * settings.dt * 1e3;
                locs(pks<=0) = [];
                pks(pks<=0) = [];

                plot(locs,pks,'.','MarkerSize',15);

        %         locs(pks== 1e3*(Vspike - E)) = [];
        %         pks(pks == 1e3*(Vspike - E)) =[];
                tmpks = pks;
                d = diff(locs); % finding the time interval between peaks
                Y = discretize(d,bin); % splitting the time intervals into the appropriate bins
                try
                    pks(1)= [];
                catch
                end
                occ = mode(Y);
                curVal = nan(sum(Y==occ),length(bin));
                for j = 1 : nanmax(Y)
                    if sum(Y==j) ~=0
                        curVal(1:sum(Y==j),j) = pks(Y==j);
                    end
                end
                val = [val ;curVal];
                pks = tmpks;

            end

        end
%         V = neuronGroup{i}.Post_Synaptic_Potential;
%         plot((0:length(V(1,:))-1)*settings.dt*1e3,V*1e3,'.');
        xlabel('Time [msec]');
        ylabel('EPSP [mV]');
        title('EPSP as a function of time');
        set(gca,'XLim',[0 settings.tmax*1e3])
        set(gca,'YLim',[0 (params.Vth - params.E)*1e3]);

        saveas(h,fullfile(outPath,'EPSP as a function of time.bmp'));
        a = sum(~isnan(val),1) < 5;
        val(:,a) = nan;
        h = figure('Visible','Off');
        hold all;
        Y = nanmean(val);
        Y(isnan(Y)) = [];
        X = (0:length(Y)-1)*binWidth;
        
        Er = nanstd(val)/(sqrt(sum(~isnan(val(:)))));
        Er(isnan(Er)) = [];
        bar(X,Y)
        errorbar(X,Y,Er,'.');

%         X = (0:length(nansum(val))-1)*binWidth;
%         Y = nansum(val);
%         bar(X,Y);

        xlabel('PSP interval width [msec]');
        ylabel('PSP [mV]');
        title('Average PSP as a function of interval width');
        set(gca,'XLim',[-binWidth/2 10*binWidth]);
        set(gca,'YLim',[0 (params.Vth - params.E)*1e3]);
        saveas(h,fullfile(outPath,'PSP vs interval.bmp'));


%         h = figure('Visible','Off');
%         imagesc(~V);
        tVec = 0:settings.dt:settings.tmax-settings.dt;
        h = plotRaster(logical(neuronGroup{i}.Post_Synaptic_Spike_Times),tVec*1e3);
        colormap(gray);
        xlabel('Time');
        ylabel('Neuron');
        title('Post Synaptic Raster Plot');
        saveas(h,fullfile(outPath,'Post Synaptic Raster Plot.bmp'));
        handler = h;

        stepSize = settings.SW/settings.dt;
        for n = 1 : (numel(tVec) - stepSize)
            spikes = neuronGroup{i}.Post_Synaptic_Spike_Times(:,n:n+stepSize-1);
            Sync_Indicator(i) = Sync_Indicator(i) + sum(spikes(:));
        end
%         Sync_Indicator(i) = Sync_Indicator(i)/params.num_of_postsynaptic_neurons;
%         Sync_Indicator(i) = nanmean(SyncIdx)/params.num_of_postsynaptic_neurons; % normalizing the synchronicity indicator by the number of neurons
    %     disp(['The synchronicity index is: ' num2str(nanmean(SyncIdx))]);
        Sync_Indicator(i) = Sync_Indicator(i)/sum(neuronGroup{i}.Post_Synaptic_Spike_Times(:));
        if settings.one_on_one
            u = neuronGroup{i}.u(1,:);
            x = neuronGroup{i}.x(1,:);
            h = figure('Visible','Off');
            subplot(311);
            V = neuronGroup{i}.Post_Synaptic_Potential;
            plot((0:length(V(1,:))-1)*settings.dt*1e3,V(1,:)*1e3);
            ylabel('Post synaptic potential');
            subplot(312);
            hold on;
            tVec = (0:length(V(1,:))-1)*settings.dt*1e3;
            plot(tVec(u~=0),u(u~=0),'-o');
            plot(tVec(x~=0),x(x~=0),'-o');
            set(gca,'XLim',[0 settings.tmax*1e3]);
            legend('u','x','Location','East');
            ylabel('Synaptic Variables');

            subplot(313);
            V = neuronGroup{i}.Pre_Synaptic_Spikes(1,:);
            tVec = find(V)*settings.dt*1e3;
            stem(tVec,V(V~=0),'.');
            ylabel('Pre synaptic spikes');
            set(gca,'XLim',[0 settings.tmax*1e3]);

            xlabel('Time [msec]');
            saveas(h,fullfile(outPath,'u_x_Potential.bmp'));
        end
        if settings.one_on_one_multi
            u = neuronGroup{i}.u(idx(1),:);
            x = neuronGroup{i}.x(idx(1),:);
            h = figure('Visible','Off');
            subplot(311);
            V = neuronGroup{i}.Post_Synaptic_Potential;
            plot((0:length(V(1,:))-1)*settings.dt*1e3,V(1,:)*1e3);
            ylabel('Post synaptic potential');
            subplot(312);
            hold on;
            tVec = (0:length(V(1,:))-1)*settings.dt*1e3;
            plot(tVec(u~=0),u(u~=0),'-o');
            plot(tVec(x~=0),x(x~=0),'-o');
            set(gca,'XLim',[0 settings.tmax*1e3]);
            legend('u','x','Location','East');
            ylabel('Synaptic Variables');
            
            plotRaster(logical(neuronGroup{i}.Post_Synaptic_Spike_Times),tVec,h);

            xlabel('Time [msec]');
            saveas(h,fullfile(outPath,'u_x_Potential.bmp'));
        end
        
        
        S = neuronGroup{i}.Connectivity_Ratio;
        h = visualise_connectivity(S);
        saveas(h,fullfile(outPath, 'Connectivity.bmp'));
        
        S_AC = neuronGroup{i}.AC_Connectivity_Ratio;
        h = visualise_connectivity(S_AC);
        saveas(h,fullfile(outPath,'AC Connectivity.bmp'));
        close all;
        
        
        
        
        
    end
    h = figure;
    bar(1:numel(Sync_Indicator),Sync_Indicator);
    axis tight;
    if numel(params.pAC) > 1
        xlabel('Probability of AC connections');
        ylabel('Synchronicity indicator');
        title('Synchronicity indicator as a function of connection probability');
        set(gca,'XTickLabel',num2cell(params.pAC));
    else
        xlabel('Model Mode');
        ModelModes = {'Strong Depression','Depression','Facilitation-Depression','Facilitation','Strong Facilitation'};
        set(gca,'XTickLabel',ModelModes,'XTickLabelRotation',45);
        title('Synchronicity indicator as a function of model mode');
        ylabel('Synchronicity indicator');
    end
    saveas(h,'Figures\Sync Indicator.bmp');
    
    
    h = figure('Visible','Off');
    num_of_post_synaptic_spikes = num_of_post_synaptic_spikes/params.num_of_postsynaptic_neurons;
    bar(1:numel(num_of_post_synaptic_spikes),num_of_post_synaptic_spikes);
    axis tight;
    if numel(params.pAC) > 1
        xlabel('Probability of AC connections');
        ylabel('Number of post synaptic spikes/Number of post synaptic neurons');
        title('Number of post synaptic spikes as a function of connection probability');
        set(gca,'XTickLabel',num2cell(params.pAC));
    else
        xlabel('Model Mode');
        ModelModes = {'Strong Depression','Depression','Facilitation-Depression','Facilitation','Strong Facilitation'};
        set(gca,'XTickLabel',ModelModes,'XTickLabelRotation',45);
        title('Number of post synaptic spikes as a function of model mode');
        ylabel('Number of post synaptic spikes/Number of post synaptic neurons');
        
    end
    saveas(h,fullfile(pwd,'Figures','Number of post synaptic spikes.bmp'));
    
    
    h = figure('Visible','Off');
    sub_treshold_EPSP = sub_treshold_EPSP/params.num_of_postsynaptic_neurons;
    bar(1:numel(sub_treshold_EPSP),sub_treshold_EPSP);
    axis tight;
    if numel(params.pAC) > 1
        xlabel('Probability of AC connections');
        ylabel('Total amount of sub threshold EPSPs/Number of post synaptic neurons');
        title('Total amount of sub threshold EPSPs as a function of connection probability');
        set(gca,'XTickLabel',num2cell(params.pAC));
    else
        xlabel('Model Mode');
        ModelModes = {'Strong Depression','Depression','Facilitation-Depression','Facilitation','Strong Facilitation'};
        set(gca,'XTickLabel',ModelModes,'XTickLabelRotation',45);
        title('Total amount of sub threshold EPSPs as a function of model mode');
        ylabel('Total amount of sub threshold EPSPs/Number of post synaptic neurons');
        
    end
    saveas(h,fullfile(pwd,'Figures','Total amount of sub threshold EPSPs.bmp'));
    
    
    

        
end
