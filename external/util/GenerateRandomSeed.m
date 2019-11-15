function rndSeed = GenerateRandomSeed(varargin)
% function GenerateRandomSeed([nextra])
%
% Returns a random, machine specific integer to use as a random seed. 
if (~isempty(varargin))
  nextra = varargin{1};
  if (nextra<=0)
    error('SetRandomSeed: Argument must be > 0.');
  end
else
  nextra = 1;
end
[ret, sys] = system('hostname');
clk = clock;
clk = round(clk(end)*1000);
rndSeed = SimpleIntegerHash(mod(sum(sys),256) + mod(clk,1024) + nextra);
