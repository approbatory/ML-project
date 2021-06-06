t_ = tic;
max_rate = 0.25;
g_field = @(x, m, s) max_rate.*exp(-((x-m)./s).^2 ./ 2);

n_cells = 50e3;
track_len = 120;
field_std = 3;

x = 0:0.1:track_len; % cm
m = sort(rand(n_cells, 1) * track_len);

R = g_field(x, m, field_std);

rate_func = @(x) g_field(x, m, field_std);
response_sample = @(x) poissrnd(rate_func(x));

t = seconds(seconds(0):seconds(1/20):hours(0.5)); %s
max_v = 30; %cm/s
x_traj = track_len .* sin(max_v./track_len .* t).^2;

noise_sigma = 10;
input_noise = noise_sigma*randn(size(x_traj));

X_noise = response_sample(x_traj + input_noise);

tracesEvents.simulated_noise = X_noise.';
tracesEvents.position = [x_traj(:), 0*x_traj(:)];
savefile = 'Mouse-2099-20210606_000000-linear-track-TracesAndEvents.mat';
save(savefile, 'tracesEvents');

d_sim_noise = DecodeTensor({savefile, 'Mouse2099'}, 'simulated_noise');

%%

d_neu = 10;
n_neu = unique([1, d_neu:d_neu:n_cells, n_cells]);

progressbar('neu...');
for n_ix = 1:numel(n_neu)
    [~, mse(n_ix)] = d_sim_noise.basic_decode(false, n_neu(n_ix), [], my_algs('lda'));
    [~, mse_sh(n_ix)] = d_sim_noise.basic_decode(true, n_neu(n_ix), [], my_algs('lda'));
    progressbar(n_ix/numel(n_neu));
end

figure;
plot(n_neu, 1./mse); hold on;
plot(n_neu, 1./mse_sh);
xlabel 'Num neurons'
ylabel '1/MSE (cm^{-2})'
toc(t_);