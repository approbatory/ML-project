function [r,a,b,c,d] = ComputeRandIndex(L1, L2)
% [r,a,b,c,d] = ComputeRandIndex(L1, L2)
%
% Given two sets of cluster labels L1 and L2 for the same set of
% points, computes the rand index measuring their similarity, and
% returns it in R.  The index is 0 if the clusterings in L1 and L2
% disagree completely, and 1 if they agree completely, and
% intermediate otherwise.
%
% The function also returns the constituents of the rand index computation:
%
% A: The number of pairs that were in the same cluster for both clusterings.
% B: The number of pairs that were in different clusters for both clusterings.
% C: The number of pairs that were in the same cluster for L1, but different in L2.
% D: The number of pairs that were in the same cluster for L2, but different in L1.
%
% R is then computed as
%
% R = (A+B)/(A+B+C+D)

if (~isvector(L1) | ~isvector(L2))
  error('Inputs must be vectors.');
end

L1 = L1(:);
L2 = L2(:);

if (numel(L1) ~= numel(L2))
  error('Label vectors must be the same size.');
end

e1 = 1 - pdist([L1 L1],'hamming');
e2 = 1 - pdist([L2 L2],'hamming');

% Gotta use this instead of 'dot' since dot doesn't like logicals,
% and the '~'s below create logicals.
f = @(x,y) sum(x.*y);

a = f(e1,e2);
b = f(~e1,~e2);
c = f(e1,~e2);
d = f(~e1,e2);

r = (a+b)/(a+b+c+d);


