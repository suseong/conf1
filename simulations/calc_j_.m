function acc = calc_j_(a,t)

x = zeros(size(t));

order = 13;
chad = zeros(order);
for k=1:order-3
    chad(3+k,k) = k*(k+1)*(k+2);
end

for i=1:size(t,2)
    t_ = calc_t(t(i),order);
    acc(i) = a'*chad*t_';
end

end