function bitsRecovered = decMpsk(M, MpskSymbols, symbolMap)
    % Pre-allocate recovered symbol indices
    numSymbols = length(MpskSymbols);
    recoveredIndices = zeros(numSymbols, 1);

    % Map each received symbol to the nearest constellation point
    for k = 1:numSymbols
        [~, idx] = min(abs(MpskSymbols(k) - symbolMap)); % nearest constellation point
        recoveredIndices(k) = idx; 
    end

    % Convert Gray code indices back to binary indices
    if M == 2
        % For BPSK, simply map +1 -> 1, -1 -> 0
        bitsRecovered = (real(MpskSymbols) > 0);
    else
        % Reverse Gray coding
        grayIndices = recoveredIndices - 1;          % MATLAB 1-index adjustment
        binIndices = zeros(size(grayIndices));
        for n = 1:length(grayIndices)
            g = grayIndices(n);
            b = 0;
            mask = g;
            while mask > 0
                b = bitxor(b, mask);
                mask = bitshift(mask, -1);
            end
            binIndices(n) = b;
        end

        % Convert integer symbol indices to bits
        bitsPerSymbol = log2(M);
        bitsRecovered = de2bi(binIndices, bitsPerSymbol, 'left-msb'); 
        bitsRecovered = reshape(bitsRecovered.', [], 1);  % column vector
    end
end