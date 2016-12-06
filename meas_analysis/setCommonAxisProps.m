function setCommonAxisProps()

    alw = 0.75;    % AxesLineWidth
    fsz = 10;      % Fontsize
    lw = 0.75;%1.25;      % LineWidth
    msz = 6.0;       % MarkerSize
    
%    grid on
    set(gca,'XGrid','on')
    set(gca,'XMinorGrid','off')
    set(gca,'YGrid','on')
    set(gca,'YMinorGrid','off')
    set(gca,'GridAlpha',0.25)
    set(gca,'MinorGridAlpha',0.4)
    set(gca,'Fontsize',fsz)
    set(gca,'LineWidth',alw);
    set(gca,'FontName','TimesRoman')
    
    % set the line properties
    hline = get(gca,'Children');
    for h = hline(:)'
        if strcmp(h.Type,'line')
            h.LineWidth = lw;
            h.MarkerSize = msz;
        end
    end
end

