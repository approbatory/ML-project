%function spikest = mlspike_nocal(trace, dt)
%if ~exist('dt', 'var')
%    dt = 1/20;
%end
function spikest = mlspike_nocal(trace, dt, varargin)
p = inputParser;
p.addRequired('trace', @isvector);
p.addOptional('dt', 1/20, @isscalar);
p.addParameter('show', false, @islogical);
p.parse(trace, dt, varargin{:});


% parameters
% (get the default parameters)
par = tps_mlspikes('par');
% (indicate the frame duration of the data)
par.dt = dt;
% (set physiological parameters)
par.a = 0.07; % DF/F for one spike
par.tau = 1; % decay time constant (second)
par.saturation = 0.1; % OGB dye saturation
par.hill = 1.7;
% (set noise parameters)
par.finetune.sigma = 0.125;%.02; % a priori level of noise (if par.finetune.sigma
                          % is left empty, MLspike has a low-level routine 
                          % to try estimating it from the data
par.drift.parameter = .015; % if par.drift parameter is not set, the 
                           % algorithm assumes that the baseline remains
                           % flat; it is also possible to tell the
                           % algorithm the value of the baseline by setting
                           % par.F0
% (do not display graph summary)
par.dographsummary = false;

% spike estimation
calcium = trace + 1;
[spikest, fit, drift] = spk_est(calcium,par);


if p.Results.show
    figure;
    spk_display(dt,{[], spikest},{calcium fit drift});
end
end