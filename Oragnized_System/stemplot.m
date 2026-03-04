function fig_n = stemplot(fig_n, x, y, ttl, xlbl, ylbl, varargin)
    fig_n = fig_n + 1;
    figure(fig_n);
    stem(x,y,varargin{:});
    title(ttl);
    xlabel(xlbl); ylabel(ylbl);
    grid on;
end