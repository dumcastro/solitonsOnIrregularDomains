function [sigma] = lmb2sig(lambda) %gets stnd deviation so that gaussian pulse has lambda as eff wavelength

threshold = 0.01;
sigma = sqrt(-lambda^2/(8*log(threshold)));

end