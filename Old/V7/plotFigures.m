function plotFigures (neuronGroup,params,settings)
    Sync_Indicator = zeros(numel(params.pAC),1);
    for i = 1 : numel(params.pAC)
        mkdir(['Figures\pAC ' num2str(params.pAC(i))]);
        
        V = neuronGroup{i}.Post_Synaptic_Potential;
        I = neuronGroup{i}.Post_Synaptic_Current;
        h = figure('Visible','Off');
        suptitle('Post Synaptic Potential')
        idx = find((mean(V,2) - params.E) > 1e-4, 4);
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

        

        saveas(h,'Figures\Post Synaptic Potential.bmp');

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
        saveas(h,'Figures\Synaptic Current.bmp');
        
        
        V = neuronGroup{i}.Post_Synaptic_Potential - params.E;
        
        
        
        
        
        V(V == params.Vspike-params.E) = nan;

        h = figure('Visible','Off');
        hold all;
        binWidth = 50;
        bin = 0:binWidth:settings.tmax*1e3;
        val = nan(1,length(bin));
        for m = 1 : size(V,1)
            [pks,locs] = findpeaks(V(m,:)*1e3); % This will locate the location and size of the EPSPs.
            if ~isempty(pks)
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
        hold off;
%         V = neuronGroup{i}.Post_Synaptic_Potential;
%         plot((0:length(V(1,:))-1)*settings.dt*1e3,V*1e3,'.');
        xlabel('Time [msec]');
        ylabel('EPSP [mV]');
        title('EPSP as a function of time');

        saveas(h,['Figures\pAC ' num2str(params.pAC(i)) '\EPSP as a function of time.bmp']);
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

        axis tight
        xlabel('PSP interval width [msec]');
        ylabel('PSP [mV]');
        title('Average PSP as a function of interval width');
        axis tight;
        saveas(h,['Figures\pAC ' num2str(params.pAC(i)) '\PSP vs interval.bmp']);


%         h = figure('Visible','Off');
        V = neuronGroup{i}.Post_Synaptic_Potential - params.E;
        V(V<params.Vspike-params.E) = 0;
        V(V == params.Vspike-params.E) = 1;
%         imagesc(~V);
        tVec = 0:settings.dt:settings.tmax-settings.dt;
        h = plotRaster(logical(V),tVec*1000/2);
        colormap(gray);
        xlabel('Time');
        ylabel('Neuron');
        title('Post Synaptic Raster Plot');
        saveas(h,['Figures\pAC ' num2str(params.pAC(i)) '\Post Synaptic Raster Plot.bmp']);

        
        for n = 1 : (numel(tVec) - settings.SW/settings.dt)
            spikes = neuronGroup{i}.Spike_Times(n:n+settings.SW/settings.dt);
            Sync_Indicator(i) = Sync_Indicator(i) + sum(spikes(:));
        end
        Sync_Indicator = Sync_Indicator/params.num_of_postsynaptic_neurons;
%         Sync_Indicator(i) = nanmean(SyncIdx)/params.num_of_postsynaptic_neurons; % normalizing the synchronicity indicator by the number of neurons
    %     disp(['The synchronicity index is: ' num2str(nanmean(SyncIdx))]);
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
            saveas(h,'Figures\u_x_Potential.bmp');
        end
        
        
        S = neuronGroup{i}.Connectivity_Ratio;
        h = visualise_connectivity(S);
        saveas(h,['Figures\pAC ' num2str(params.pAC(i)) '\Connectivity.bmp']);
        
        S_AC = neuronGroup{i}.AC_Connectivity_Ratio;
        h = visualise_connectivity(S_AC);
        saveas(h,['Figures\pAC ' num2str(params.pAC(i)) '\AC Connectivity.bmp']);
        close all;
        
        
        
        
    end
    h = figure;
    bar(params.pAC,Sync_Indicator);
    axis tight;
    xlabel('Probability of AC connections');
    ylabel('Synchronicity indicator');
    title('Synchronicity indicator as a function of connection probability');
    saveas(h,'Figures\Sync Indicator.bmp');
end
