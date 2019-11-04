function Y = TocToMatrix(X,d)
% Y = TocToMatrix(X, d)
%
% Converts the Toc Matrix with column dimensions d into a 4-D matrix.
%
% Y(t,o,c,:) is the data for the specified trial, odor and cell.

Y = reshape(X, [], d(1), d(2), d(3));
Y = permute(Y, [2 3 4 1]);