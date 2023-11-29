function [value, Freq] = histRate(x)
%   Copyright xiezhh，2010.3.8

x = sort(x(:));
x1 = diff(x);
x1(end + 1) = 1;
x1 = find(x1);
value = x(x1);
x1 = [0; x1];
Freq = diff(x1);