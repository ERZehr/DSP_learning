function ROM = genPolyphaseROM(h, L, direction)

    % Compute number columns in the ROM matrix
    row_size = length(h) / L;

    % Reshape into L polyphase rows
    base = reshape(h, L, row_size);

    % Direction selection
    if strcmp(direction, 'counter-clockwise')
        ROM = base;
    elseif strcmp(direction, 'clockwise')
        ROM = flipud(base);
    else
        error('direction must be ''clockwise'' or ''counter-clockwise''');
    end

end
