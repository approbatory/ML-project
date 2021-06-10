function s = fastpnb

s.name = 'Fast Poisson NB';
s.pre = @(x)x;
s.train = @pnb_train;
s.test = @pnb_test;
