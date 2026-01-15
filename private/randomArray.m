function rand_array = randomArray(n, min_val, max_val)
    rand_array = min_val + (max_val - min_val) * rand(1, n);
end