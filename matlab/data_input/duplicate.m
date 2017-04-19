% duplicate
%
% Given a labelled dataset, this will return a new dataset with each class
% repeated according to a given vector of desired proportions. For example,
% duplicate(X, Y, 0.1, [1, 1, 2] will return the same dataset with data points
% corresponding to the 3rd class doubled (with noise). duplicate(X, Y, 0.1, [2,
% 2, 2] will double the entire dataset and add noise to the added data points.
% 0s in the proportions vector are ignored; there will always be at least one
% copy of the original classes.
%
% PARAMETERS:
%   X -- A N x K matrix of N data points each of dimension K, arranged in rows
%
%   Y -- A N x 1 vector of class labels, where Y(i) gives the class label of
%        data point X(i, :)
%
%   noise_var -- The variance of the Gaussian noise that will be added to the
%                duplicated data points.
%
%   proportions -- A 1 x C vector, where there are C class labels, where
%                  proportions(i) gives the relative proportions of the
%                  original class i that is desired in the new dataset.
%
% RETURNS
%   Xd, Yd -- The new X and Y for the duplicated data
%
% Written by Jesse Hellemn at Rice University, 2016
function [Xd, Yd] = duplicate(X, Y, noise_var, proportions)
    Xd = X;
    Yd = Y;

    % Duplicate every class
    for c = 1:size(proportions, 2)
        c_idx = Yd == c;
        Xc = Xd(c_idx, :);
        Yc = Yd(c_idx);
        
        % Start at 2 because 1 copy already exists
        for dup = 2:proportions(c)
            Xd = [Xd ; Xc + noise_var*randn(size(Xc))];
            Yd = [Yd ; Yc];
        end
    end
end
