function plot_daysets_end_decoding(daysets, res, varargin)
p = inputParser;
p.addRequired('daysets', @iscell);
p.addRequired('res', @isstruct);
%p.addParameter('matching', false, @islogical);
p.addParameter('save_to', '', @ischar);
p.addParameter('suppress', false, @islogical);
p.parse(daysets, res, varargin{:});

matching = res.was_matched;
filling = res.filling;
mode = res.mode;
res = res.res;

if p.Results.suppress
    set(0,'DefaultFigureVisible','off');
end

colors = 'xbxgxrx';
circle_begin = 0.4374;
for m_ix = 1:numel(daysets)
    figure;
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
        title([mode ' decoder error on matching cells (on ' filling ')'], 'Interpreter', 'none');
    else
        title([mode ' decoder error on all available cells (on ' filling ')'], 'Interpreter', 'none');
    end
    xlabel('fraction of turn completed');
    ylabel('Linear SVM decoding error');
    ylim([0 0.5]);
    line([circle_begin circle_begin], ylim, 'Color', 'k',...
        'DisplayName', 'entering center');
    if ~isempty(p.Results.save_to)
        if ~exist(p.Results.save_to, 'dir')
            mkdir(p.Results.save_to);
        end
        if matching
            fname = sprintf('%s_%s_%s_matched', daysets{m_ix}(1).day(1:5), mode, filling);
        else
            fname = sprintf('%s_%s_%s_unmatched', daysets{m_ix}(1).day(1:5), mode, filling);
        end
        print(fullfile(p.Results.save_to, fname), '-dpng');
    end
end
if p.Results.suppress
    set(0,'DefaultFigureVisible','on');
end
end