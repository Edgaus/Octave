function bowing = bowing_function(x, b0, x0, n, s)
    bowing = b0 ./ ( 1 + (x ./ x0).^n ).^s;
end