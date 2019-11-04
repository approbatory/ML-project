function Q = ComputeSubplotPositions(numRows, numCols, Prc, hMarginLeft, hMarginRight, hSpacing, vMarginTop, vMarginBottom, vSpacing)
% Q = ComputeSubplotPositions(numRows, numCols, Prc, 
%                             hMarginLeft, hMarginRight,  hSpacing, 
%                             vMarginTop,  vMarginBottom, vSpacing)
%
% Computes the positions of the subplots whose unit positions are
% specified in Prc, taking into account the margin and spacing
% values provided. If Prc is empty, the positions for 
% numRows x numCols unit plots are computed.
%
% The subplot positions are first computed by breaking the figure
% into a [numRows x numPlot] grid of unit subplots, taking the
% spacing and margin information provided into account. Subplot
% positions are specified by indicating overlap with these
% unitplots. The subplots are specified as a cell array of
% cells. A cell specifying a subplot is a pair of vectors
% specifying the row and column indices of the unit plots it should
% overlap. For example, {{1,1}, {1,2:3}, {2,1}, {3,1}, {2:3,2}} 
% specifies 5 subpplots as below
%
%                             1 2 2
%                             3 5 5
%                             4 5 5
%
% Q is a [nSubPlots x 4] matrix specifying the positions of the
% subplot as [left bottom width height] vectors suitable for
% passing to the subplot function. To see the positions of the
% subplots, call ShowComputedSubplotPositions.
%
% Examples:
%
% Q = ComputeSubplotPositions(1,1,{{1,1}},0.01,0.01,0.01,0.02,0.02,0.02);
%
% Computes the positions for a single subplot with horizontal
% margins of 0.01 and vertical margins of 0.02.
%
% Q = ComputeSubplotPositions(2,2,{{1,1},{1,2},{2,1},{2,2}},0.01,0.01,0.02,0.01,0.01,0.02)
%
% Compute the positions of four subplots, placed at the top left,
% top right, bottom left and bottom right of the figure. Horizontal
% and vertical margins are both 0.01, horizontal spacing is 0.02,
% vertical margins are 0.01, vertical spacing is 0.02.
%
% Q = ComputeSubplotPositions(2,3,{{1,1:2},{1,3},{2,1},{2,2:3}},0.01,0.01,0.02,0.01,0.01,0.02)
%
% Computes the position of four subplots, the first a [1 x 2] plot
% at the top left corner, the second a [1 x 1] plot at the top
% right corner, the third a [1 x 1] plot at the bottom left corner
% and the fourth a [1 x 2] plot at the bottom right corner.
%
% Q = ComputeSubplotPositions(7,1,0.01,0.01,0.02,0.01,0.01,0.02)
%
% Computes the positions for a stack of 7 equally sized plots.


% First do some range checking to make sure the positioning
% arguments are OK.

if (numRows<1)
  error('Number of rows must be >= 1.');
end

if (numCols<1)
  error('Number of columns must be >= 1.');
end

occupied = zeros(numRows, numCols);

if (isempty(Prc))
  k = 1;
  for i = 1:numRows
    for j = 1:numCols
      Prc{k} = {i,j};
      k = k+1;
    end
  end
end

for i = 1:numel(Prc)
  rows = Prc{i}{1};
  cols = Prc{i}{2};
  
  if (min(rows)<1)
    error('Subplot %d: found at least one row index < 1.', i);
  end
  
  if (max(rows)>numRows)
    error('Subplot %d: found at least one row index > max.', i);
  end
  
  if (min(cols)<1)
    error('Subplot %d: found at least one column index < 1.', i);
  end
  
  if (max(cols)>numCols)
    error('Subplot %d: found at least one column index > max.', i);
  end
  
  if (numel(rows)>1)
    if (unique(diff(rows))~=1)
      error('Subplot %d: rows in range must be adjacent.',i);
    end
  end

  if (numel(cols)>1)
    if (unique(diff(cols))~=1)
      error('Subplot %d: columns in range must be adjacent.',i);
    end
  end
    
  currentOccupants = occupied(rows,cols);
  currentOccupants = unique(currentOccupants(:));
  nzCurrentOccupants = currentOccupants(find(currentOccupants));
  nzCurrentOccupants = nzCurrentOccupants(:)';
  if (~isempty(nzCurrentOccupants))
    error('Position of subplot %d conflicts with subplot(s) %s',i,num2str(nzCurrentOccupants));
  end
  
  occupied(rows,cols) = i;
end
    

% Figure out the unit sizes.

unitHeight = (1-vMarginTop - vMarginBottom-vSpacing*(numRows-1))/numRows;
unitWidth  = (1-hMarginLeft- hMarginRight -hSpacing*(numCols-1))/numCols;

unitLeft   = ((1:numCols)-1)*(unitWidth+hSpacing) + hMarginLeft;
unitBottom = ((numRows:-1:1)-1)*(unitHeight+vSpacing)+ vMarginBottom;


% Assemble the subplot positions from the units.

Q = zeros(numel(Prc),4);
for i = 1:numel(Prc)
  rmin = min(Prc{i}{1});
  rmax = max(Prc{i}{1});
  cmin = min(Prc{i}{2});
  cmax = max(Prc{i}{2});
  
  Q(i,1) = unitLeft(cmin);
  Q(i,2) = unitBottom(rmax);
  Q(i,3) = unitWidth *(cmax-cmin+1)+hSpacing*(cmax-cmin);
  Q(i,4) = unitHeight*(rmax-rmin+1)+vSpacing*(rmax-rmin);
end

