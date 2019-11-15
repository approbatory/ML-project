function visualize_resids(X,Xest,md,trial_map)

% calculate error and resid
resid = Xest-X;
err = resid.^2;

% break down errors and resids across trial types
Rsq.east = [];
Rsq.west = [];
Rsq.correct = [];
Rsq.incorrect = [];
for k = 1:size(trial_map,1)
    ek = err(:,:,k);
    xk = X(:,:,k);
    rsq = 1 - sum(ek(:))/sum((xk(:)-mean(xk(:))).^2);
    
    d = trial_map(k,1);
    kd = trial_map(k,2);
    
    start = md.day(d).trials(kd).start;
    if md.day(d).trials(kd).correct
        correctness = 'correct';
    else
        correctness = 'incorrect';
    end
    Rsq.(start) = [Rsq.(start); rsq];
    Rsq.(correctness) = [Rsq.(correctness); rsq];
end

% make the plot
figure()

subplot(2,2,1)
histfit(resid(:))
xlabel('residual')
legend('observed residual','gaussian fit')
title('histogram of all residuals')

subplot(2,2,3); hold on
C = [1,2,4,5];
mn_rsq = [mean(Rsq.east), mean(Rsq.west), mean(Rsq.correct), mean(Rsq.incorrect)];
bar(C,mn_rsq,'r')
title('Model fit east/west starts, correct/incorrect trials (bars show means)')
ylabel('R^2')

N = [length(Rsq.west), length(Rsq.east), length(Rsq.correct), length(Rsq.incorrect)];

i = 1;
for rsq = {Rsq.west, Rsq.east, Rsq.correct, Rsq.incorrect}
    n = N(i);
    c = C(i);
    plot(c*ones(n,1)+0.1*randn(n,1),rsq{1},'.k')
    i = i+1;
end
set(gca,'xtick',[1,2,4,5],'xticklabels',{'west','east','correct','error'})

subplot(1,2,2)
E = log(squeeze(sum(err,2)));
%[~,Es1] = sort(sum(E,1),'descend');
[~,Es2] = sort(sum(E,2),'descend');
image(E(Es2,:),'CDataMapping','scaled')
xlabel('trial number')
ylabel('neurons')
title('heatmap of log sum-of-squared residuals for each neuron on each trial')
