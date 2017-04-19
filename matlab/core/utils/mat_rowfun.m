% mat_rowfun
%
% Same thing as Matlab's built in rowfun.m, but works on matrices instead of
% tables. This function will apply a function to every row of a matrix; the
% given function must take in a row vector and return a row vector. The
% dimension of the returned vector does not need to be the same size as the
% input vector, but it must be consistent across all rows of the matrix.
%
% This function is actually implemented with a for loop; there may be more
% efficient implementations possible. The dimensions of the output matrix are
% assumed to be the output dimension of func applied to the first row of A.
%
% Arguments:
%
% A                     The matrix to call functions on.
%
% func                  The function to call on each row of A. This must take
%                       in a row vector as its only parameter and output a row
%                       vector as its first return value.
%
% Written by Jesse Hellemn, Spring 2017 at Rice University
function B = mat_rowfun(A, func)

    % First detect the dimensionality of the resulting matrix
    first_row = func(A(1,:));

    % Enforce that the result is a row
    if numel(first_row) ~= max(size(first_row))
        error('func must produce a row vector when applied to another row vector')
    end

    % Allocate memory
    B = zeros(size(A, 1), numel(first_row));

    % Put first row in
    B(1,:) = first_row;

    % Call the function on every row
    for i = 2:size(A, 1)
        B(i,:) = func(A(i,:));
    end

end
