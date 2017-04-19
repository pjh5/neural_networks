% confusion_matrix
%
% Returns the confusion matrix (as a matrix) of classifying desired as output.
%
% PARAMETERS:
% 
% desired   The true class labels, in 1-in-C encoding, arranged in rows.
% output    The predicted class labels, in 1-in-C encoding, arranged in rows.
%   
% Written by Jesse Hellemn at Rice University, 2016
function cmat = confusion_matrix(desired, output)
    [~,prediction] = classification_accuracy(desired, output);
    cmat = desired' * double(prediction);
end
