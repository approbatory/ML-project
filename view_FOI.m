function view_FOI(frame_of_interest, pos)
%VIEW_FOI Summary of this function goes here
%   Detailed explanation goes here
figure
for i = 1:size(frame_of_interest,1)
    c = pos{i};
    plot(c(:,1), c(:,2), 'b.');
    hold on
    xlim([0 1]);
    ylim([0 1]);
end
colors = {'r', 'g', 'c', 'y', 'm', 'k'};
for i = 1:size(frame_of_interest,1)
    c = pos{i};
    hold on
    xlim([0 1]);
    ylim([0 1]);
    for j = 1:size(frame_of_interest,2)
        plot(c(frame_of_interest(i,j),1), c(frame_of_interest(i,j),2), [colors{mod(j-1,length(colors))+1} 'x']);
    end
end
end

