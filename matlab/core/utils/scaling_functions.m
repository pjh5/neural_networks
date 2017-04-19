% scaling_functions
%
% Returns two functions, mapping from the data space of Y to the unit
% square, and the inverse mapping.
%
% ARGUMENTS:
%   Y - N x K numeric array of N data points, each with K dimensions. Every
%       row is a separate data point.
%
%   scale -- The data will be mapped onto [-scale, scale]. Whether each
%            dimension will be scaled individually is up to the next
%            parameter. This should really just be 1 in most cases.
%
%   by_column -- Whether to scale each column individually or not. If the
%                columns are independent, then they should really be scaled
%                independently (this should be true).
%
% RETURNS:
%   f_scale -- A function that, when given any matrix with K columns, will
%              scale it according to the parameters above. Generally, scaling
%               
%
%   f_unscale -- A function that, when given any matrix with K columns, will
%                perform the inverse mapping from the scaled space back to the
%                original data space. This is exactly the inverse of f_scale.
%                Be careful not to unscale data twice, this makes no checks to
%                see if the input matrix lives in the scaled space.
%
% Written by Jesse Hellemn at Rice University, Fall 2016
function [f_scale, f_un] = scaling_functions(Y, scale, by_column)

    % Vectorized functions for column scaling
    if (by_column)
        rangeY = range(Y,1);
        midY = (max(Y,[],1) + min(Y,[],1)) / 2;

        % Don't divide by range if a column is only one value
        rangeY(rangeY == 0) = 1;

        % Return scaling functions
        f_scale = @(y) 2*scale * ...
            bsxfun(@rdivide, ...
                   bsxfun(@minus, y, midY), ...
                   rangeY);
        f_un = @(y) bsxfun(@plus, ...
                           bsxfun(@times, y, rangeY) / (2*scale), ...
                           midY);

    % Non-Vectorized functions for global scaling
    else
       rangeY = max(max(Y)) - min(min(Y));
       midY = (max(max(Y)) + min(min(Y))) / 2;

       % Don't scale by range if range is 0
       if rangeY == 0
           rangeY = 1;
       end

       % Return scaling functions
       f_scale = @(y) 2*scale * (y - midY) / rangeY;
       f_un = @(y) y * rangeY / (2*scale) + midY;
    end
end
