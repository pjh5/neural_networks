% read_fcps
%
% Use this function to read in one of the FCPS datasets.
% FCPS stands for the Fundamental Clustering Problems Suite, first created for
%    Ultsch, A.: Clustering with SOM: U*C, In Proc. Workshop on Self-Organizing
%    Maps, Paris, France, (2005), pp. 75-82
% The datasets are free for download from the Philipps University website at:
% https://www.uni-marburg.de/fb12/arbeitsgruppe
%
%
% ARGUMENTS:
%   dataset_name -- The name of the dataset, as a string. For example "Atom"
%                   or "Lsun". Note that all of the datasets begin with a
%                   capitalized letter.
%
%   proportions -- (Optional) A vector of size C x 1, where the dataset
%                  specified by the first parameter has C classes in it. If
%                  present, proportions will be passed to duplicate() to
%                  duplicate the classes of the dataset . See the comment in
%                  duplicate.m for more detail.
%
%   noise_var -- (Optional, but required if the 2nd parameter is present) The
%                variance of the gaussian noise to add to duplicated data
%                points. See the comments in duplicate.m for more detail.
%
% RETURNS:
%   X -- An N x K numeric matrix of N data points, each of dimension K. This
%        is the raw input data, and does not include class labels.
%   Y -- An N x 1 numeric array, where the ith element is the class label of
%        the ith datapoint (the ith row of X)
%
% Written by Jesse Hellemn at Rice University, 2016
function [X, Y] = read_fcps(dataset_name, proportions, noise_var)

    % Declare path to the FCPS directory, where all datasets can be found
    base_path = strcat('~/pjh5/neural_nets/data/FCPS/01FCPSdata/', dataset_name);

    % Data is stored in .lrn files
    % The first column is just an index, which is redundant
    X = load(strcat(base_path, '.lrn'));
    X = X(:, 2:end);

    % Class labels are stored in .cls files
    % The first column is just an index, which is redundant
    Y = load(strcat(base_path, '.cls'));
    Y = Y(:, 2:end);

    % If duplicate requested, then return duplicated data
    if nargin > 1
        [X, Y] = duplicate(X, Y, noise_var, proportions);
    end
end
