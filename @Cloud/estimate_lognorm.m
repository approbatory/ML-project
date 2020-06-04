sigma = linspace(0, 5, 500);


D = zeros(size(sigma));
D_sem = D;

D_unbiased = D;
D_unbiased_sem = D;

N = 100;
samp = 1e2;
dens = @(x) (sum(x).^2 ./ sum(x.^2)) ./ size(x,1);
dens_unbiased = @(x) exp(-var(log(x)));
progressbar('sigma...');
for i = 1:numel(sigma)
    X = exp(sigma(i)*randn(N, samp));
    temp = dens(X);
    temp_unbiased = dens_unbiased(X);
    
    D(i) = mean(temp);
    D_sem(i) = sem(temp);
    
    D_unbiased(i) = mean(temp_unbiased);
    D_unbiased_sem(i) = sem(temp_unbiased);
    progressbar(i/numel(sigma));
end

%%
%errorbar(sigma, D, D_sem);
figure;
plot(sigma, D);
hold on;
plot(sigma, D_unbiased);
xlabel '\sigma'
ylabel 'Log normal signal density'

hold on;

plot(sigma, exp(-sigma.^2));

%plot(sigma, exp(-sigma.^2)*(1-1/N) + 1/N);
legend 'Simulated' 'Unbiased' 'Model exp(-\sigma^2)' %'Model2'

title(sprintf('N = %d',N));

%%
figure;
semilogy(sigma, D .* exp(sigma.^2));
hold on;
plot(sigma, exp(sigma.^2));