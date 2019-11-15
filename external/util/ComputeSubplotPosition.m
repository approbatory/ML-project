function Q = ComputeSubplotPosition(leftMargin, rightMargin, topMargin, bottomMargin)
% Q = ComputeSubplotPosition(leftMargin, rightMargin, topMargin, bottomMargin)
%
% A wrapper for ComputeSubplotPositions for the case of a single subplot.

Q = ComputeSubplotPositions(1,1,[],leftMargin,rightMargin,0, topMargin, bottomMargin, 0);
