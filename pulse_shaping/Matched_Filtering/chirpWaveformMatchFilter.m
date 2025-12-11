% Pulse Wave
close all; clear; clc;

% Generate Waveform
[pulseWaveform, t] = waveformGenerator('chirp', 6e6, 5e-5, 500e3);
figure(1)
stem(t, real(pulseWaveform));
title('Pulse');
xlabel('Time (s)'); ylabel('Magnitude');


% apply pulse shaping
[shapedWaveform, returnVars] = pulseShape('raisedcosine', pulseWaveform, [10, 2, 5, 0.2]);
figure(2)
stem(real(shapedWaveform));
title('Shaped Pulse Waveform');
xlabel('Time (s)'); ylabel('Magnitude');


% Generate Gauss noise vector
noiseLength = 50000;
gaussNoise = randn(noiseLength,1)';
figure(3)
plot(gaussNoise);


%generate a random index in the gauss noise for my signal to "come back"
startDelayLength = 30000;
startIndex = randi(startDelayLength);


% generate the signal that comes back
% recievedSignal = gauss noise + signal at the random time that it comes back
delayedSignal = [zeros(1,startIndex) shapedWaveform zeros(1, noiseLength - (startIndex + length(shapedWaveform)))];
idx = find(delayedSignal ~= 0, 1, 'first'); % check where the signal starts
fprintf('Actual Delay Index = %d\n', idx);
% Plot the delayed signal
figure(4);
plot(delayedSignal);
title('Delayed Signal');
xlabel('Index'); ylabel('Magnitude');
grid on;


% Calculate the recieved data
% recievedSignal = gauss noise + signal at the random time that it comes back
recievedSignal = gaussNoise + delayedSignal;
% Plot the delayed function
figure(5);
plot(recievedSignal);
title('Recieved Signal');
xlabel('Index'); ylabel('Magnitude');
grid on;


% Calculate the autocorrelation of the sequence
autocorreleation = zeros(1,length(recievedSignal)-length(shapedWaveform) + 1);
for k = 1:length(autocorreleation)
    recievedChunk = recievedSignal(k:k+length(shapedWaveform)-1);
    autocorreleation(k) = sum(shapedWaveform .* recievedChunk);
end

% Find the max of autocorrelation
[autocorreleationMaxValue, autocorreleationMaxIndex]  = max(autocorreleation);
fprintf("Calculated Delay Index = %d", autocorreleationMaxIndex);

% Plot the autocorrelation function
figure(6);
plot(autocorreleation);
title('Autocorrelation of sent and recieved data');
xlabel('Index'); ylabel('Normalized Magnitude');
grid on;