function traj = visualize_trial_trajectories(cpd, meta)
% VISUALIZE_TRIAL_TRAJECTORIES(cpd, trial_idx)

% number of trials to plot
K = length(meta.start)

% number of factors in the model
nf = size(cpd.factors.trial,2);

% enumerate all combinations of factors
allcombs = combnk(1:nf,3);
ncomb = nchoosek(nf,3);

% calculate explanatory power of all factor combinations
pwr = zeros(ncomb,1);
for c = 1:ncomb
    comb = allcombs(c,:);
    pwr(c) = sum(abs(cpd.lambda(comb)));
end

% sort factor combs by power (most interesting plots are shown first)
[pwr,cmborder] = sort(pwr,'descend')

% iterate over within-trial factor combinations
figure()
for c = 1:ncomb
    clf()

    % current combination of within-trial factors
    comb = allcombs(cmborder(c),:);
    c1 = comb(1);
    c2 = comb(2);
    c3 = comb(3);

    % plot all trials
    for k = 1:K
        % trial index
        Bk = cpd.factors.time * diag(cpd.factors.trial(k,:));

        x = Bk(:,c1);
        y = Bk(:,c2);
        z = Bk(:,c3);

        % color start location
        subplot(2,2,1)
        hold on
        switch meta.start{k}
            case 'east' % blue
                plot3(x, y, z, '-', 'linewidth', 0.5, 'color', [0.0 0.7 1.0])
            case 'west' % red
                plot3(x, y, z, '-r', 'linewidth', 0.5)
            otherwise % black
                plot3(x, y, z, '-k', 'linewidth', 0.5)
        end

        % color end location
        subplot(2,2,2)
        hold on
        switch meta.end{k}
            case 'north' % blue
                plot3(x, y, z, '-', 'linewidth', 0.5, 'color', [0.0 0.7 1.0])
            case 'south' % red
                plot3(x, y, z, '-r', 'linewidth', 0.5)
            otherwise % black
                plot3(x, y, z, '-k', 'linewidth', 0.5)
        end

        % color correct/incorrect
        subplot(2,2,3)
        hold on
        switch meta.correct{k}
            case '1' % blue
                plot3(x, y, z, '-', 'linewidth', 0.5, 'color', [0.0 0.7 1.0])
            case '0' % red
                plot3(x, y, z, '-r', 'linewidth', 0.5)
            otherwise % black
                plot3(x, y, z, '-k', 'linewidth', 0.5)
        end

        % color day/session number
        subplot(2,2,4)
        hold on
        switch meta.day{k}
            case '12' % blue
                plot3(x, y, z, '-', 'linewidth', 0.5, 'color', [0.0 0.7 1.0])
            case '13' % red
                plot3(x, y, z, '-r', 'linewidth', 0.5)
            otherwise % black
                plot3(x, y, z, '-k', 'linewidth', 0.5)
        end

    end

    subplot(2,2,1)
    title(sprintf('Factors (%i,%i,%i) start location',c1,c2,c3))
    grid on
    
    subplot(2,2,2)
    title(sprintf('Factors (%i,%i,%i) end location',c1,c2,c3))
    grid on

    subplot(2,2,3)
    title(sprintf('Factors (%i,%i,%i) correct/incorrect',c1,c2,c3))
    grid on
    
    subplot(2,2,4)
    title(sprintf('Factors (%i,%i,%i) days',c1,c2,c3))
    grid on

    pause
end
