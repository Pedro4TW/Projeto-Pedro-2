margin = 0.5;
Xlim = [ min( cell2mat({country_bounds.X}) )-margin;max( cell2mat({country_bounds.X}) )+margin ];
Ylim = [ min( cell2mat({country_bounds.Y}) )-margin;max( cell2mat({country_bounds.Y}) )+margin ];
set( gca,'XTick',[],'YTick',[],'XLim',Xlim,'YLim',Ylim,'Box','on' )