function session_filterer(n, d)
n_sess = 239;
my_section = ceil((1:n_sess)./n_sess.*d) == n;
my_inds = find(my_section);

addpath decoding/
progressbar('filtering...');
fname = sprintf('session_filterer_%d_%d.txt', n, d);
fid = fopen(fname,'w');
%m_e = zeros(1, 239);
counter = 0;
for i = my_inds
    try
        d = DecodeTensor.cons(i);
        m_e = d.basic_decode(false, [], []);
        fprintf('%d %f\n', i, m_e);
        fprintf(fid, '%d %f\n', i, m_e);
    catch e
        fprintf('%d 1000\n', i);
        fprintf(fid, '%d 1000\n', i);
    end
    counter = counter + 1;
    progressbar(counter/numel(my_inds));
end