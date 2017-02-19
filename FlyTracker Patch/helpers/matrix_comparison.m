
% Determine if the values of personal and environmental features have not
% been changed by the transformations of the patch implemented for
% calculating the relative features of 3+ flies.
% 
% This is achieved by taking a copy of the dataset before relative features
% have been computed, and compare them with the same columns of the new,
% healed file.
% 
% Using the "isequaln" method allows comparing each individual value in
% these matrices, treating nans and missing values as comparable and equal.
% This method will return a single boolean value of 1 (True) if, and only if,
% absolutely ALL values are identical in the two matrices; otherwise,
% returns 0 (False).

% Load the matrices to compare.
before = feat_before.data(:,:,:);
after = feat_after.data(:,:,1:9);

% Check the size of each matrix to make sure they are the same size across
% all dimensions (matrices must be exactly the same size in order to be
% compared by the "isequaln" method).
size_before = size(before)
size_after = size(after)

% Actually compare the individual values in each of the matrices.
is_equal = isequaln(before, after)