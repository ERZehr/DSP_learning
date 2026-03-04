function firpm_h = gen_firpm_first_order_h(pass_type, wp, ws, Rp, Rs)
    dev_pass = (10^(Rp/20)-1)/(10^(Rp/20)+1);
    dev_stop = 10^(-Rs/20);
    switch pass_type
        case 'low'
            f = [wp/pi ws/pi];
            a = [1 0];
            dev = [dev_pass dev_stop];
        case 'high'
            f = [ws/pi wp/pi];
            a = [0 1];
            dev = [dev_stop dev_pass];
        otherwise
            f = [wp/pi ws/pi];
            a = [1 0];
            dev = [dev_pass dev_stop];
    end
    [N,fo,ao,w] = firpmord(f,a,dev);
    isEven = mod(N,2) == 0;
    N = N + isEven; % round to nearest odd N
    firpm_h = firpm(N,fo,ao,w);
end