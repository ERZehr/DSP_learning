function ROM = barkerROM(h,L)
    row_size = length(h)/L;
    ROM = reshape(h, L, row_size).';
end