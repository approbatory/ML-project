function spikest = mlspike_yescal(trace, dt, varargin)
p = inputParser;
p.addRequired('trace', @isvector);
p.addOptional('dt', 1/20, @isscalar);
p.addParameter('show', false, @islogical);
p.parse(trace, dt, varargin{:});

calcium = trace + 1;


amin = .04;
amax = .1;
taumin = .4;
taumax = 1.6;
%sigmamin = .005;
%sigmamax = .05;

% parameters
% (get default parameters and set frame duration)
pax = spk_autocalibration('par');
pax.dt = dt;
% (set limits for A and tau)
pax.amin = amin;
pax.amax = amax;
pax.taumin = taumin;
pax.taumax = taumax;
% (set saturation parameter)
pax.saturation = 0.1;
% (give real values to the algorithm for display purposes - they obviously won't be used in the estimation!)
%pax.realspikes = spikes;
%pax.reala = a;
%pax.realtau = tau;
% (when running MLspike from spk_autocalibratio, do not display graph summary)
pax.mlspikepar.dographsummary = false;
pax.display = 'none';

% perform auto-calibration
[tauest, aest, sigmaest] = spk_autocalibration(calcium,pax);
fprintf('Estimated tau=%f\ta=%f\tsigma=%f\n', tauest, aest, sigmaest);

% parameters
par = tps_mlspikes('par');
par.dt = dt;
% (use autocalibrated parameters)
par.a = aest;
par.tau = tauest;
par.finetune.sigma = sigmaest;
% (the OGB saturation and drift parameters are fixed)
par.saturation = 0.1;
par.hill = 1.7;
par.drift.parameter = .015/5;
% (do not display graph summary)
par.dographsummary = false;

% spike estimation
[spikest, fit, drift] = spk_est(calcium,par);

if p.Results.show
    figure;
    spk_display(dt,{[], spikest},{calcium fit drift});
end
end