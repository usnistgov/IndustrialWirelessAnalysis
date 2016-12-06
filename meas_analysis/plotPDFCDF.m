function h = plotPDFCDF(X, N, metric, units)
        X = X(~isnan(X));
        [X_counts,X_centers] = hist(X,N);
        X_probs = cumsum(X_counts/sum(X_counts));     
        X_pdf = [0 diff(X_probs)];
        
        h = figure();  
        
        % PDF 
        yyaxis left, harea = area(X_centers, X_pdf, 'FaceAlpha',0.5);
        str = 'Pr. $$ \hat{X} $$'; 
        ylabel(str,'Interpreter','Latex');   
        harea.EdgeColor = 'none';
        
        % CDF
        yyaxis right, plot(X_centers,X_probs)
        ylim([0 1])        
        str = 'Pr. $$ \hat{X} < {X} $$';
        ylabel(str,'Interpreter','Latex'); 
        
        % xlabel
        xlstr = sprintf('%s, $${X}$$ (%s)', metric, units);
        xlabel(xlstr,'Interpreter','Latex');
        
        drawnow
end