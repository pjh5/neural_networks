% range
%
% Returns the range max(A) - min(A) in the given dimension. If no dimension is
% given, the default is 1.
function [result] = range(A, dim)

    % Optional dim parameter
    if nargin < 2
        dim = 1
    end

    result = max(A, [], dim) - min(A, [], dim);

end
