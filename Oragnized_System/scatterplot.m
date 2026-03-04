function fig_n = scatterplot(fig_n, x, y, ttl, xlbl, ylbl, varargin)
    fig_n = fig_n + 1;
    figure(fig_n);
    scatter(x,y,varargin{:});
    title(ttl);
    xlabel(xlbl); ylabel(ylbl);
    grid on;
end