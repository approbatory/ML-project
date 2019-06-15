name = 'effect_of_mindist_on_curves_2022';
figure;
dbfile = 'decoding_with_mindist.db';
conn = sqlite(dbfile);

%carry out sql query using session_id to avoid future ambiguity
% sessions with mindist tests: DecodeTensor.qid([1 2 4 5 6 7])
mindist_set = 0:5:15;
m_i_indices = [1 2 4 5 6 7];
res = db_knife(conn, 'neuron_series', DecodeTensor.qid(m_i_indices),...
    'settings_set', {'unshuffled', 'shuffled'},...
    'mindist_set', mindist_set);
%%
my_lines = [];
for i = m_i_i:m_i_i
    for j = 1:2
        for k = [1 3 4]
            setts = {'unshuffled', 'shuffled'};
            my_lines = [my_lines DecodingPlotGenerator.errors_plotter(res.n{i,j,k}, res.m{i,j,k}, res.e{i,j,k}, setts{j}, 'index', k, 'DisplayName', sprintf('Dist. = %d px', mindist_set(k)))];
        end
    end
end
xlabel 'Number of cells'
ylabel '1/MSE (cm^{-2})'
xlim([0 500]);
text(300, 0.25/6, 'Unshuffled', 'Color', 'blue');
text(400, 0.6/6, 'Shuffled', 'Color', 'red');
legend(my_lines, 'Location', 'northwest');
%at the end
conn.close;
d_ = DecodeTensor(m_i_indices(m_i_i));
title(d_.mouse_name);
return;
%%
%todo: ball and stick plots showing no difference in curve params after
%fits

%changes to N?
figure;
subplot(2,2,1);
ballnstick('Unrestricted', '10 px minimal distance', res.N_fit(:, 1, 1), res.N_fit(:, 1, 3), res.N_fit_errb(:, 1, 1), res.N_fit_errb(:, 1, 3));
set(gca, 'YScale', 'log');
ylabel('N fit (neurons)');
title('Unshuffled');

%changes to I_0?
subplot(2,2,2);
ballnstick('Unrestricted', '10 px minimal distance', 1e4*res.I_0_fit(:, 1, 1), 1e4*res.I_0_fit(:, 1, 3), 1e4*res.I_0_fit_errb(:, 1, 1), 1e4*res.I_0_fit_errb(:, 1, 3));
%set(gca, 'YScale', 'log');
ylabel('I_0 fit (m^{-2}/neuron)');
title('Unshuffled');

%changes to N?
subplot(2,2,3);
ballnstick('Unrestricted', '10 px minimal distance', res.N_fit(:, 2, 1), res.N_fit(:, 2, 3), res.N_fit_errb(:, 2, 1), res.N_fit_errb(:, 2, 3));
set(gca, 'YScale', 'log');
ylabel('N fit (neurons)');
title('Shuffled');

%changes to I_0?
subplot(2,2,4);
ballnstick('Unrestricted', '10 px minimal distance', 1e4*res.I_0_fit(:, 2, 1), 1e4*res.I_0_fit(:, 2, 3), 1e4*res.I_0_fit_errb(:, 2, 1), 1e4*res.I_0_fit_errb(:, 2, 3));
%set(gca, 'YScale', 'log');
ylabel('I_0 fit (m^{-2}/neuron)');
title('Shuffled');

%function res = db_knife(conn, series_type, mice_session_set, settings_set, mindist_set, aggregate_flag, diff_flag)
function res = db_knife(varargin)
p = inputParser;
p.addRequired('conn', @(x)true);
p.addRequired('series_type', @ischar);
p.addRequired('mice_session_set', @(x) (ischar(x) || iscell(x)));
p.addParameter('settings_set', 'unshuffled', @(x) (ischar(x) || iscell(x)));
p.addParameter('mindist_set', 0, @isnumeric);
p.addParameter('aggregate_flag', false, @islogical);
p.addParameter('diff_flag', false, @islogical);

p.parse(varargin{:});
r = p.Results;

if strcmp(r.series_type, 'neuron_series')
    primary_field = 'NumNeurons';
    secondary_field = 'DataSize';
elseif strcmp(r.series_type, 'datasize_series')
    primary_field = 'DataSize';
    secondary_field = 'NumNeurons';
else
    error('Only ''neuron_series'' or ''datasize_series'' is allowed as series_type');
end

basic_fetcher = @(sess, sett, mindist)...
    sprintf('select %s, MSE from decoding where SessionID = ''%s'' and Setting = ''%s'' and MinDist = %f and %s = (select max(%s) from decoding where SessionID = ''%s'' and Setting = ''%s'' and MinDist = %f);',...
    primary_field, sess, sett, mindist, secondary_field, secondary_field, sess, sett, mindist);

if ischar(r.mice_session_set)
    sessions = {r.mice_session_set};
else
    sessions = r.mice_session_set;
end

if ischar(r.settings_set)
    settings = {settings_set};
else
    settings = r.settings_set;
end

for i_s = 1:numel(sessions)
    my_sess = sessions{i_s};
    for j_set = 1:numel(settings)
        my_setting = settings{j_set};
        for k_md = 1:numel(r.mindist_set)
            my_mindist = r.mindist_set(k_md);
            [res.n{i_s, j_set, k_md}, res.m{i_s, j_set, k_md}, res.e{i_s, j_set, k_md},...
                res.N_fit(i_s, j_set, k_md), res.I_0_fit(i_s, j_set, k_md),...
                res.N_fit_errb(i_s, j_set, k_md), res.I_0_fit_errb(i_s, j_set, k_md)] =...
                aggregate_one_curve(my_sess, my_setting, my_mindist);
        end
        %correcting fit bias, to use the same support for n values
        n_upper_limit = min(cellfun(@max, res.n(i_s, j_set, :)));
        for k_md = 1:numel(r.mindist_set)
            my_mindist = r.mindist_set(k_md);
            [~, ~, ~,...
                res.N_fit(i_s, j_set, k_md), res.I_0_fit(i_s, j_set, k_md),...
                res.N_fit_errb(i_s, j_set, k_md), res.I_0_fit_errb(i_s, j_set, k_md)] =...
                aggregate_one_curve(my_sess, my_setting, my_mindist, n_upper_limit);
        end
    end
end

    function [n, m, e, N_fit, I_0_fit, N_fit_errb, I_0_fit_errb] = aggregate_one_curve(sess, sett, mindist, n_upper_limit)%TODO fill out
        if ~exist('n_upper_limit', 'var')
            n_upper_limit = Inf;
        else
            fprintf('Using upper limit for n : %d\n', n_upper_limit);
        end
        db_ret = r.conn.fetch(basic_fetcher(sess, sett, mindist));
        [n_rows, n_cols] = size(db_ret);
        n_vals = cell2mat(db_ret(:,1));
        imse_vals = 1./cell2mat(db_ret(:,2));
        n = unique(n_vals);
        for i = 1:numel(n)
            my_n = n(i);
            my_samples = imse_vals(n_vals == my_n);
            m(i) = mean(my_samples);
            e(i) = std(my_samples)./sqrt(length(my_samples))*norminv((1+0.95)/2);
        end
        
        filt_n = double(n(n <= n_upper_limit));
        filt_m = m(n <= n_upper_limit);
        [fitresult, gof] = createFit_infoSaturation(filt_n(:), filt_m(:));
        N_fit = fitresult.N;
        I_0_fit = fitresult.I_0;
        conf_errbs = diff(confint(fitresult))/2;
        I_0_fit_errb = conf_errbs(1);
        N_fit_errb = conf_errbs(2);
    end

    function [n_agg, m_agg, e_agg] = aggregate_many_curves(n_cell, m_cell)%TODO fill out
    end

end

%%knifelang:
% neuron_series [2022 2023 2024] {unshuffled shuffled}
% neuron_series [2022 2023 2024 2029] {unshuffled diagonal} -agg -diff
% neuron_series [2022] {unshuffled shuffled} <0 5 10 15 20 25>