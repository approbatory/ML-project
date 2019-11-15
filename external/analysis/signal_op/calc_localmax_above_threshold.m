function localmax = calc_localmax_above_threshold(threshold,data)
%function calc_localmax_above_threshold(threshold,data)
%
%find the local maxima in a sequence of points provided that the maxima are 
%above a specified threshold
%
%Hakan Inan (Jan 15)
%
data_binary = data>threshold; %binarize
%vector to contain local maxima
localmax = [];
k=1;
while k<length(data_binary)
    if data_binary(k)==1 % start a sequence of 1's
        begin_seq = k; 
        while (data_binary(k)==1)&&(k<length(data_binary)), k=k+1; end
        end_seq = k-1;
        %after the end of sequence calculate max and add it to localmax
        seq = data(begin_seq:end_seq); %magnitudes of the events in the sequence
        [maxval,dum] = max(seq);
        idx_max = dum+begin_seq-1;
        if maxval>threshold*1.1 % reject very small peaks
            localmax = [localmax,idx_max]; %#ok
        end        
    else %skip through 0's
        k = k+1;
    end
end