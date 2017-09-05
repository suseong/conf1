function acc = calc_a_(a,t)

x = zeros(size(t));

order = 13;
chad = zeros(order);
for k=1:order-2
    chad(2+k,k) = k*(k+1);
end

for i=1:size(t,2)
    t_ = calc_t(t(i),order);
    acc(i) = a'*chad*t_';
end

end