function [h] = genPulseFilter(rolloff, span, L, type)
    h = rcosdesign(rolloff, span, L, type);
    h = h(1:end-1); % trim one coeff so we have a rectangular ROM 
    h = h / max(h);
end