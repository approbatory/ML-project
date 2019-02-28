n = 0:5:500;
I_s = n;
I_us = n./(1 + n/300);
I_d = n./(1 + n/100);
figure; plot(n, I_s, 'r'); hold on; plot(n, I_us, 'b'); plot(n, I_d, 'm');
xlabel 'Number of cells'
ylabel 'Place information'
plot(n, n.*0 + 300, 'b:')
plot(n, n.*0 + 100, 'm:')
legend Shuffled Unshuffled 'Correlation-insensitive'
legend Shuffled Unshuffled 'Correlation-insensitive' Location NorthWest
xlim([0 500]);
%axis equal
set(gca, 'YTick', [])
legend off
figure_format