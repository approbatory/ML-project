paths = ...
{'/home/omer/_distalopto_traces/per-nvoke/m834/m834-0710';
'/home/omer/_distalopto_traces/per-nvoke/m834/m834-0712';
'/home/omer/_distalopto_traces/per-nvoke/m834/m834-0711';
'/home/omer/_distalopto_traces/per-nvoke/m833/m833-0704';
'/home/omer/_distalopto_traces/per-nvoke/m833/m833-0703';
'/home/omer/_distalopto_traces/per-nvoke/m833/m833-0705';
'/home/omer/_distalopto_traces/mcherry/m893/m893-1028';
'/home/omer/_distalopto_traces/mcherry/m893/m893-1027';
'/home/omer/_distalopto_traces/mcherry/m893/m893-1026';
'/home/omer/_distalopto_traces/mcherry/m895/m895-1107';
'/home/omer/_distalopto_traces/mcherry/m895/m895-1106';
'/home/omer/_distalopto_traces/mcherry/m895/m895-1108';
'/home/omer/_distalopto_traces/mcherry/m891-120trials/m891-1024';
'/home/omer/_distalopto_traces/mcherry/m891-120trials/m891-1023';
'/home/omer/_distalopto_traces/mcherry/m891-120trials/m891-1025';
'/home/omer/_distalopto_traces/mcherry/m892/m892-1016';
'/home/omer/_distalopto_traces/mcherry/m892/m892-1018';
'/home/omer/_distalopto_traces/mcherry/m892/m892-1017';
'/home/omer/_distalopto_traces/mcherry/m894/m894-1104';
'/home/omer/_distalopto_traces/mcherry/m894/m894-1102';
'/home/omer/_distalopto_traces/mcherry/m894/m894-1103';
'/home/omer/_distalopto_traces/icpp/m911/m911-1209';
'/home/omer/_distalopto_traces/icpp/m911/m911-1207';
'/home/omer/_distalopto_traces/icpp/m911/m911-1208';
'/home/omer/_distalopto_traces/icpp/m911/m911-1206';
'/home/omer/_distalopto_traces/icpp/m902/m902-1210';
'/home/omer/_distalopto_traces/icpp/m902/m902-1211';
'/home/omer/_distalopto_traces/icpp/m902/m902-1209';
'/home/omer/_distalopto_traces/enphr/m953/m953-0123';
'/home/omer/_distalopto_traces/enphr/m953/m953-0124';
'/home/omer/_distalopto_traces/enphr/m953/m953-0122';
'/home/omer/_distalopto_traces/enphr/m954/m954-0122';
'/home/omer/_distalopto_traces/enphr/m954/m954-0124';
'/home/omer/_distalopto_traces/enphr/m954/m954-0123';
'/home/omer/_distalopto_traces/enphr/m955/m955-0121';
'/home/omer/_distalopto_traces/enphr/m955/m955-0122';
'/home/omer/_distalopto_traces/enphr/m955/m955-0123';
'/home/omer/_distalopto_traces/enphr/m952/m952-0124';
'/home/omer/_distalopto_traces/enphr/m952/m952-0122';
'/home/omer/_distalopto_traces/enphr/m952/m952-0123';
'/home/omer/_distalopto_traces/enphr/m932/m932-0124';
'/home/omer/_distalopto_traces/enphr/m932/m932-0123';
'/home/omer/_distalopto_traces/enphr/m932/m932-0125';
'/home/omer/_distalopto_traces/enphr/_sham/m952-0125';
'/home/omer/_distalopto_traces/enphr/_sham/m932-0126';
'/home/omer/_distalopto_traces/enphr/_sham/m953-0125';};


for i = 1:6
    [n_inhib(i), n_disinhib(i), n_total(i)] = process_opto_session(paths{i});
    fprintf('%s\tinhib:%d\tdisibhib:%d\n', paths{i}, n_inhib(i), n_disinhib(i));
end

%%
for i = 1:6
    if n_inhib(i) < 0
        fprintf('%s\terror code %d\n', paths{i}, n_inhib(i));
    else
        fprintf('%s\tinhib: %d\tdisinhib: %d\n', paths{i}, n_inhib(i), n_disinhib(i));
    end
    pause
end