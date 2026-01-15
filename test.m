x = linspace(-10*pi,10,1000);

f = exp(-x.^2);

K = 42;

f_shift = circshift(f,K);

