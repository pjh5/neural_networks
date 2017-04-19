% read_wdbc
%
% Returns the X and Y for the wdbc dataset
function [X_wdbc, Y_wdbc] = read_wdbc()
    formatspec = ...
            '%f%C%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
    wdbc = readtable(...
            '~/pjh5/neural_nets/data/breast_cancer_wisconsin/wdbc.csv', ...
            'Delimiter', ',', 'Format', formatspec);
    X_wdbc = table2array(wdbc(:, 3:end));
    Y_wdbc = 1 + (table2array(wdbc(:, 2)) == 'B');
end
