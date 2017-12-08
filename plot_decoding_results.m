function plot_decoding_results(decoding_results, errtype, errdesc)
if nargin == 1
    errtype = 'dist';
    errdesc = 'bin distance';
end
figure;
for ix = 1:length(decoding_results)
    subplot(1,3,ix);
    label = decoding_results(ix).label;
    res = decoding_results(ix).res;
    algs = cell(1,length(res));
    div_types = cell(1,3);
    for i = 1:length(res)
        algs{i} = res(i).alg.name;
        div = res(i).division;
        for j = 1:length(div)
            div_types{j} = div(j).desc;
            out = div(j).out;
            for k = 1:length(out)
                if strcmp(out(k).err_type, errtype)
                    d(j,i) = out(k).test_err;
                end
            end
        end
    end
    bar(categorical(div_types), d);
    title(label);
    ylabel(errdesc);
    legend(algs);
end

end