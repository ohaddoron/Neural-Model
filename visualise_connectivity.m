function h = visualise_connectivity(S)

    h = figure('Visible','Off');
    Ns = size(S,1); % Number of pre-synaptic cells 
    Nt = size(S,2); % Number of post-synaptic cells
    hold on;
    plot(zeros(Ns,1),1:Ns,'LineWidth',5);
    plot(ones(Nt,1),1:Nt,'LineWidth',5);
    for i = 1 : Ns
        a = find(S(i,:));
        for j = 1 : length(a)
            line([0,1],[i,a(j)],'LineWidth',0.0001,'Color','k');
        end
    end
    axis off;
    
end