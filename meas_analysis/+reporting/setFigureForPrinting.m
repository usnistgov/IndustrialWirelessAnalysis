function setFigureForPrinting(h)
    width=3.5; height=3.5;
    %set(gcf,'InvertHardcopy','on');
    set(h,'PaperUnits', 'inches');
    papersize = get(h, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(h,'PaperPosition', myfiguresize);
end