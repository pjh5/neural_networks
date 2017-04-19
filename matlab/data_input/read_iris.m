% read_iris
%
% Returns the X and Y (in 1-in-C encoding) of the iris dataset
function [X, Yc] = read_iris()

    % Read in the iris dataset
    iris = readtable('~/pjh5/neural_nets/data/iris.csv');

    % Remove unnecessary index column
    iris = iris(:, 2:end);

    % Convert species name to 1-in-C encoding
    iris.Species = categorical(iris.Species);
    iris.is_setosa     = iris.Species == 'setosa';
    iris.is_versicolor = iris.Species == 'versicolor';
    iris.is_virginica  = iris.Species == 'virginica';
    iris.Species = [];

    % Convert the table object to an array
    iris = table2array(iris);

    % Split into input and output
    X = iris(:, 1:4);
    Yc = iris(:, 5:7);
end
