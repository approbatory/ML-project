function view_FOI(frame_of_interest, pos)
%VIEW_FOI Plot where the frames of interest are
if isstruct(frame_of_interest)
    [frame_of_interest, pos] = find_frame_at_pos(frame_of_interest, pos);
end

%pts = linspace(0,1,200);
%[X,Y] = meshgrid(pts,pts);

%bins = bin_f(X,Y);
%disp(unique(bins(:))');
figure
%sizes = 10*mod(bins,2)+1;
%scatter(X(:),Y(:),sizes(:),bins(:));
%hold on;
for i = 1:size(frame_of_interest,1)
    c = pos{i};
    %plot(c(:,1), c(:,2), 'k.');
    scatter(c(:,1), c(:,2), 1, bin_space(c));
    hold on
    xlim([0 1]);
    ylim([0 1]);
end
colors = {'r', 'g', 'b', 'c', 'y', 'm'};
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