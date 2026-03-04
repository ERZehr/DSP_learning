function signal_out = apply_channel_noise(signal_in,gauss_noise_level,salt_and_pepper_noise_level,linear_channel_attenuation)
    gaussNoise = randn(length(signal_in),1)'*gauss_noise_level;

    % Generate Salt and Pepper noise vector
    prob = 0.05; % probability of impulse
    saltAndPepperNoise = zeros(1,length(signal_in));
    impulses = rand(1,length(signal_in)) < prob;
    saltAndPepperNoise(impulses) = (2*randi([0,1],1,sum(impulses))-1)'*salt_and_pepper_noise_level;

    total_noise = gaussNoise + saltAndPepperNoise;
    signal_out = signal_in / linear_channel_attenuation + total_noise;
end