function [I_recovered,Q_recovered] = IQ_unmixer(signal_in,carrier_frequency,sample_rate)
    % Split I and Q
    len_cos = 0:length(signal_in)-1;
    cos_lo_rec = cos(2*pi*carrier_frequency/sample_rate*len_cos);
    sin_lo_rec = sin(2*pi*carrier_frequency/sample_rate*len_cos);

    I_recovered = 2*(signal_in .* cos_lo_rec);
    Q_recovered = -2*(signal_in .* sin_lo_rec);

    % Filter at carrier frequency to eliminate undesirable trig identity components
    % we really just have to filter out 2Fs, so pass anything less then 5 and
    % cut at higher than 1. Thus, this can be a really simple PM filter
    h_IQ_recovery = gen_firpm_first_order_h('low', 0.5*pi, 0.95*pi, 1, 30);

    I_recovered = filter(h_IQ_recovery, 1, I_recovered);
    Q_recovered = filter(h_IQ_recovery, 1, Q_recovered);
end