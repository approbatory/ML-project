function [n_neu, mse, mse_sh, fr, gof] = general_decoding_curve(varargin)
p = inputParser;
p.addParameter('max_rate', 0.05, @isnumeric);
p.addParameter('model_type', "place field", @isstring);
p.addParameter('track_len', 120, @isnumeric);
p.addParameter('n_cells', 500, @isnumeric);
p.addParameter('HWHM', 15, @isnumeric);
p.addParameter('dt', seconds(1/20), @isduration);
p.addParameter('sess_time', hours(0.5), @isduration);
p.addParameter('max_v', 30, @isnumeric);
p.addParameter('noise_sigma', 10, @isnumeric);
p.addParameter('response_type', "poisson", @isstring);
p.addParameter('approx_samp', 20, @isnumeric);
p.addParameter('n_reps', 16, @isnumeric);
p.addParameter('n_min', 10, @isnumeric);
p.addParameter('progress', true, @islogical);
p.addParameter('linear_only', false, @islogical);
p.parse(varargin{:});
opt = p.Results;

max_rate = opt.max_rate;
model_type = opt.model_type;
track_len = opt.track_len;
n_cells = opt.n_cells;
HWHM = opt.HWHM;
dt = opt.dt;
sess_time = opt.sess_time;
max_v = opt.max_v;
noise_sigma = opt.noise_sigma;
response_type = opt.response_type;
approx_samp = opt.approx_samp;
n_reps = opt.n_reps;
n_min = opt.n_min;
progress = opt.progress;
linear_only = opt.linear_only;

switch model_type
    case "place field"
        field_std = 2*HWHM / 2.355; % cm, 15cm HWHM
        m = sort(rand(n_cells, 1)) * track_len;
        rate_func = @(x) max_rate.*...
            exp(-((x-m)./field_std).^2./2);
    case "rate coding"
        m = sort(2*rand(n_cells, 1)-1);
        flipped = m < 0;
        rate_func = @(x) max_rate./track_len.*...
            max(0, (x.*(1-flipped) + (track_len-x).*flipped) .* abs(m));
end

t = seconds(seconds(0):dt:sess_time); %s
x_traj = track_len .* sin(max_v./track_len .* t).^2;
position = [x_traj(:), 0*x_traj(:)];


input_noise = noise_sigma*randn(size(x_traj));

switch response_type
    case "poisson"
        response_sample = @(x) poissrnd(rate_func(x));
        alg = fastpnb;
    case "gaussian"
        output_noise = 0;%0.25;
        response_sample = @(x) rate_func(x) + output_noise.*randn(n_cells, numel(x));
        alg = fastdlda;
end

X_noise = response_sample(x_traj + input_noise);

dt_opts = DecodeTensor.default_opt;
[~,~,tr_s,tr_e,tr_dir,tr_bins,~] = DecodeTensor.new_sel(position, dt_opts);

data_tensor = DecodeTensor.construct_tensor(X_noise.', tr_bins, dt_opts.n_bins, tr_s, tr_e);


n_neu = unique(round(logspace(log10(n_min),log10(n_cells),approx_samp)));
NN = numel(n_neu);

if progress
    %progressbar('reps...', 'neu...');
    %hbar = parfor_progressbar(n_reps*NN,'Computing...');
    WaitMessage = parfor_wait(n_reps*NN, 'ReportInterval', ceil(n_reps*NN/100));
else
    WaitMessage = [];
end
[mse, mse_sh] = deal(zeros(n_reps, numel(n_neu)));

parfor r_ix = 1:n_reps
    for n_ix = 1:NN

        err_res = DecodeTensor.decode_all(data_tensor,...
            tr_dir, dt_opts.bin_width, alg, n_neu(n_ix), []);

        mse(r_ix, n_ix) = err_res.MSE.unshuffled;
        mse_sh(r_ix, n_ix) = err_res.MSE.shuffled;

        if progress
            %progressbar([],n_ix/numel(n_neu));
            %hbar.iterate(1);
            WaitMessage.Send;
        end
    end
    if progress
        %progressbar(r_ix/n_reps, []);
        %hbar.iterate(1);
    end
end
if progress
    %close(hbar);
    WaitMessage.Destroy;
end

if linear_only
    imse = 1./mse;
    
    low_bound = 0.0012;
    high_bound = 0.0462;
    
    good_vals = all(imse > low_bound,1) & all(imse < high_bound,1);
    
    if any(~good_vals)
        warning('Removed some values.');
    end
    
    n_neu = n_neu(good_vals);
    mse = mse(:, good_vals);
    mse_sh = mse_sh(:, good_vals);
    
    %load('fit_stability_analysis_workspace.mat', 'sigma_vals');
    %load('fit_stability_analysis_workspace.mat', 'rmse_vals');
    %[r_sorted, ord] = unique(rmse_vals);
    %sigma_lookup = @(r) interp1(r_sorted, sigma_vals(ord), r, 'linear', "extrap");
    %converter = @(m) sigma_lookup(sqrt(m)).^2;
    %mse = converter(mse);
    %mse_sh = converter(mse_sh);
end

if numel(n_neu) > 2    
    [fr, gof] = createFit_intercept(n_neu, mean(1./mse, 1));
else
    fr = [];
    gof = [];
end