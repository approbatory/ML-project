function ind = GetSubplotIndex(whichRow, whichCol, numRows, numCols)
% ind = GetSubplotIndex(whichRow, whichCol, numRows, numCols)
%
% Returns the index of the desired subplot in a rectangular array 
% of size numRows x numCols. 
%
% Example:
%
% Q = ComputeSubplotPositions(10,5,[],0.05,0.05,0.01,0.05,0.05,0.01);
% % plot in the subplot in the 4th column of the 7th row.
% subplotp(Q, GetSubplotIndex(7,4,10,5));
% plot(1:10);
ind = (whichRow - 1)*numCols+whichCol;