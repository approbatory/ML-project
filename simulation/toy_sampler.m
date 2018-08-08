%sample gaussian variable
N = 10;

samps = zeros(1,N);
for i = 1:N
    samps(i) = randn;
end

if ~exist('samples', 'dir')
    mkdir samples
end
save(['samples/samps' timestring '.mat'], 'samps');