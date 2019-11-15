function h = SimpleIntegerHash(n)
% h = SimpleIntegerHash(n)
% 
% A simple integer hash function, as described in 
% "Introduction to Algorithms" by Cormen et al.
A = 0.5*(sqrt(5)-1);
m = 2^10;
h = floor(m*rem(n*A,1))+1;