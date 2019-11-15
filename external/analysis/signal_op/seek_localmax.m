function x = seek_localmax(trace, x)

iter = 0;
while (iter < 100) % Hard brake to prevent long/infinite loops
    if x == 1
        delta_left = -Inf;
    else
        delta_left = trace(x-1) - trace(x);
    end

    if x == length(trace)
        delta_right = -Inf; 
    else
        delta_right = trace(x+1) - trace(x);
    end

    if (delta_left <= 0) && (delta_right <= 0)
        break;
    else
        if delta_left > delta_right
            x = x - 1;
        else
            x = x + 1;
        end
    end
    iter = iter + 1;
end
