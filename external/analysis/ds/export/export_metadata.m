function meta = export_metadata(md, trial_map)
% X = EXPORT_METADATA(md, trial_map)
%
% Exports traces of all trials specified by trial_map into
% a cell array X.

    num_trials = size(trial_map,1);

    % get moving average of turn probability on each day
%     ndays = length(md.valid_days);
%     tp = cell(ndays);
%     for di = 1:ndays
%         d = md.valid_days(di);
%         tp{di} = md.day(d).est_turn_probabilities;
%     end

    % copy selected trials into lightweight cell array
    meta.start = cell(num_trials,1);
    meta.end = cell(num_trials,1);
    meta.correct = zeros(num_trials,1);
    meta.day = zeros(num_trials,1);
    meta.turn = cell(num_trials,1);
    meta.turn_prob = zeros(num_trials,1);
    
    for k = 1:num_trials
        % day and neuron indices
        d = trial_map(k,1);
        trial = md.day(d).trials(trial_map(k,2));

        % basic metadata associated with this trial
        meta.start{k} = trial.start;
        meta.end{k} = trial.end;
        meta.correct(k) = trial.correct;
        meta.day(k) = d;
        meta.turn{k} = trial.turn;

        % turn probability estimated for this trial
%         di = md.valid_days == d;
%         meta.turn_prob(k) = tp{di}(trial_map(k,2));
    end

    % mark each trial as allo vs ego-centric
    % THK: This looks largely redundant with DaySummary.get_strategy?
    meta.strategy = cell(num_trials,1);
    e0 = NaN; % trial index of last east start
    w0 = NaN; % trial index of last west start
    for k = 1:num_trials
        pk = meta.turn_prob(k);
        if strcmp(meta.start{k},'east')
            if ~isnan(w0) && ~isnan(pk)
                pl = meta.turn_prob(w0);
                if pk > 0.99 && pl > 0.99
                    meta.strategy{k} = 'ego-right';
                elseif pk < 0.01 && pl < 0.01
                    meta.strategy{k} = 'ego-left';
                elseif pk > 0.99 && pl < 0.01
                    meta.strategy{k} = 'allo-north';
                elseif pk < 0.01 && pl > 0.99
                    meta.strategy{k} = 'allo-south';
                end
            end
            e0 = k;
        elseif strcmp(meta.start{k},'west')
            if ~isnan(e0) && ~isnan(pk)
                pl = meta.turn_prob(e0);
                if pk > 0.99 && pl > 0.99
                    meta.strategy{k} = 'ego-right';
                elseif pk < 0.01 && pl < 0.01
                    meta.strategy{k} = 'ego-left';
                elseif pk > 0.99 && pl < 0.01
                    meta.strategy{k} = 'allo-south';
                elseif pk < 0.01 && pl > 0.99
                    meta.strategy{k} = 'allo-north';
                end
            end
            w0 = k;
        else
            meta.strategy{k} = 'probe';
        end
        if isempty(meta.strategy{k})
            meta.strategy{k} = 'NA';
        end
    end
end % export_metadata