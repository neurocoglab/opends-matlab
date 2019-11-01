inter = 100;
diam_new = zeros(ceil(length(diam_avr)/inter),1);
j = 1;

for i = 1 : inter : length(diam_avr)
   ee = min(i+inter-1,length(diam_avr));
   T=diam_avr(i:ee);
   diam_new(j) = mean(T);
   j = j + 1;
end