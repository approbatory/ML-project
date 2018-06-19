function bsep(X, Ys_main, Ys_sec, max_err, legs, ylab)
if length(X)~=size(Ys_main,1) ||...
        ~isequal(size(Ys_main),size(Ys_sec))
    error(['X must have an entry for each row in Y,'...
        ' Ys_main and Ys_sec must have the same size']);
end

sem = @(x) std(x)/sqrt(length(x));
Ys_main_means = cellfun(@mean, Ys_main);
Ys_main_errbs = cellfun(sem, Ys_main);
Ys_sec_means = cellfun(@mean, Ys_sec);
Ys_sec_errbs = cellfun(sem, Ys_sec);

%my_colors = 'rgbcymk';
my_colors = {[1 0.5 0.5], [0.5 1 0.5], [0.5 0.5 1], [0.5 1 1], [1 1 0.5], [1 0.5 1], [0.5 0.5 0.5]};
color_index = @(i) my_colors{mod(i-1,length(my_colors))+1};
num_subs = length(X);
for i = 1:num_subs
    subplot(1, num_subs, i);
    b = bar(Ys_main_means(i,:));
    b.FaceColor = 'flat';
    for j = 1:size(b.CData,1)
        b.CData(j,:) = color_index(j);
    end
    hold on;
    %C = bar_center_locations(Ys_main_means(i,:));
    C = 1:length(legs{i});
    errorbar(C, Ys_main_means(i,C), Ys_main_errbs(i,C), '.k');
    plot(C, Ys_sec_means(i,C), 'o');
    errorbar(C, Ys_sec_means(i,C), Ys_sec_errbs(i,C), '.k');
    %legend(legs{i}{:});
    set(gca, 'XTickLabel', legs{i});
    set(gca, 'XTick', C);
    xtickangle(45);
    title(X{i});
    ylim([0 max_err]);
    ylabel(ylab);
end


function x = bar_center_locations(y)
[d,a] = size(y);
ds = 1:d;
as = 1:a;
x = ds' + 1/(a+1.5)*(as - (a+1)/2);
end
end