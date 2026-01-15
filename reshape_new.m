function b = reshape_new(a, rows, cols)
    % First, reshape using MATLAB's column-major ordering
    temp = reshape(a, cols, rows);
    
    % Then transpose the result to get row-major ordering
    b = temp';
end