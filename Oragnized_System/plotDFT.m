function fig_n = plotDFT(fig_n,symbols,Fs)
    dft = fft(symbols);
    dft = fftshift(dft);
    dft_axis = (-length(dft)/2:length(dft)/2-1) * (Fs/length(dft));
    contplot(fig_n, dft_axis, dft, 'DFT', 'f (Hz)', 'Magnitude');
    fig_n = fig_n + 1;
end