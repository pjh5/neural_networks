% classification_accuracy
%
% Computes the predicted class label for all outputs, and then computes the
% classification accuracy. The predicted class label from an output y^ = [y^_1,
% ..., y^_{Ky}] is argmax_i({y^_i}), or, in English, the index of the maximum
% value in the output. This is essentially assuming that the outputs are 1-in-C
% encodings of some class label, and that the given output is the
% (approximation) of the conditional probabilities of every class.
%
% PARAMETERS
% desired       The true class labels as an M x C matrix, of M inputs and C
%               possible classes. The class labels should be stored as a 1-in-C
%               encoding.
%
% output        The output of the MLP, a M x C matrix of floats. This is the
%               raw output from the MLP, for M examples, arranged in rows.
%
% RETURNS
% accuracy      The classification accuracy, as a float.
%
% prediction    The output matrix converted to 1-in-C encodings of the
%               predicted class. The predicted class for each output row is the
%               index of the maximum value in that row.
%
% Written by Jesse Hellemn at Rice University, 2016
function [accuracy, prediction] = classification_accuracy(desired, output)

    [M, C] = size(desired);
    [~, pred_class] = max(output, [], 2);
    [~, real_class] = max(desired, [], 2);
    prediction = to_one_in_c(pred_class, C);
    
    accuracy = sum(pred_class == real_class) / numel(pred_class);
end
