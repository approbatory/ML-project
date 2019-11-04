function ShowComputedSubplotPositions(Q)
% ShowComputedSubplotPositions(Q)
%
% Creates subplots in the positions specified in the matrix Q. This
% matrix is typically returned by a call to
% ComputeSubplotPositions.
%
% See also: ComputeSubplotPositions.

clf;
for i = 1:size(Q,1)
  subplot('Position',Q(i,:));
  h = text(0.5,0.5,sprintf('subplot %d',i));
  set(h,'FontSize',14,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
end