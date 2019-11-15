function Y = GetSubsetOfTocMatrix(X,d,S)
% Y = GetSubsetOfTocMatrix(X,d,S)
%
% Given the toc matrix X with toc dimensions in d, returns the subset
% of the matrix corresponding to the ranges in the 3-element cell
% array S. If any of the elements of S are empty, the data for the
% full range is returned.
%
% Example: Get the data for trials [1 3 5] for PN 143 in response to odor B
%
% >> pnSpt = LoadTocSpikeTimes('rawpns');
% >> data = GetSubsetOfTocMatrix(pnSpt, [7 44 174], {[1 3 5],5,143});
% >> data(1:5,:)
% ans =
%
%    0.0772    0.1569    0.0507
%    0.7403    0.7313    0.1439
%    0.7841    0.8277    0.2549
%    0.8377    1.3506    0.4329
%    1.0759    1.9460    1.0941

Y = TocToMatrix(X,d);

for i = 1:3
  if (isempty(S{i}))
    S{i} = 1:size(Y,i);
  end
end

Y = Y(S{1},S{2},S{3},:);
Y = MatrixToToc(Y,[numel(S{1}),numel(S{2}),numel(S{3})]);


