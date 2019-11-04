function rndSeed = SetRandomSeed1(seedOffset,varargin)
% function rndSeed = SetRandomSeed1(seedOffset) Sets the random seed
% for a single stream by Generating a random integer and adding the
% specified offset to it. Returns the seed in RNDSEED. If the optional
% argument 'override' is set to 'yes', the SEEDOFFSET is set as the
% seed directly (without having a random integer added to it first)

p = inputParser;
p.addOptional('override','no');
p.parse(varargin{:});

if (isequal(lower(p.Results.override),'yes'))
  rndSeed = seedOffset;
else
  rndSeed = GenerateRandomSeed(seedOffset);
end

stream = RandStream.create('mt19937ar','seed',rndSeed);

if (ismethod(stream, 'setGlobalStream'))
  RandStream.setGlobalStream(stream);
else
  RandStream.setDefaultStream(stream);
end
