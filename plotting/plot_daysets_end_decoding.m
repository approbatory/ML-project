function plot_daysets_end_decoding(daysets, res, varargin)
p = inputParser;
p.addRequired('daysets', @iscell);
p.addRequired('res', @isstruct);
%p.addParameter('matching', false, @islogical);
p.addParameter('save_to', '', @ischar);
p.parse(daysets, res, varargin{:});

matching = res.was_matched;
res = res.res;


figure;
colors = 'xbxgxrx';
circle_begin = 0.4374;
for m_ix = 1:numel(daysets)
    figure(m_ix);
    %loop over days
    for d_ix = 1:numel(daysets{m_ix})
        %loop over points
        if isempty(daysets{m_ix}(d_ix).changing)
            continue;
        end
        %res = daysets{m_ix}(d_ix).res;
        res_ = res{m_ix}(d_ix);
        base = res_.baseline_err;
        label = daysets{m_ix}(d_ix).day;
        errorbar(res_.points, mean(res_.test_err,2),...
            std(res_.test_err, [], 2)/sqrt(size(res_.test_err,2)),...
            colors(d_ix), 'DisplayName', label);
        l = refline(0, base);%, [colors(d_ix) '-.'], 'DisplayName', [label ' baseline']);
        l.LineStyle = '--';
        l.Color = colors(d_ix);
        l.DisplayName = [label ' baseline'];
        hold on;
    end
    legend;
    if matching
        title('End decoder error on matching cells');
    else
        title('End decoder error on all available cells');
    end
    xlabel('fraction of turn completed');
    ylabel('Linear SVM decoding error');
    ylim([0 0.5]);
    line([circle_begin circle_begin], ylim, 'Color', 'k',...
        'DisplayName', 'entering center');
    if ~isempty(p.Results.save_to)
        print(fullfile(p.Results.save_to, daysets{m_ix}(1).day(1:5)), '-dpng');
    end
end
end