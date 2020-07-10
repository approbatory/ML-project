function v = vec_ind2sub(siz, ind)
ind = ind - 1;
v = zeros(1, length(siz));

for i = 1:length(siz)
    div = floor(ind / siz(i));
    rem = mod(ind, siz(i));
    
    v(i) = rem;
    ind = div;
end

v = v + 1;