function vel = calc_v_(a,t)

x = zeros(size(t));

order = 13;
chad = zeros(order);
for k=1:order-1
    chad(1+k,k) = k;
end

for i=1:size(t,2)
    t_ = calc_t(t(i),order);
    vel(i) = a'*chad*t_';
end

end