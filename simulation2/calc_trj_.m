function x = calc_trj_(a,t)

x = zeros(size(t));
order = 13;

for i=1:size(t,2)
    t_ = calc_t(t(i),order);
   x(i) = t_*a;
end

end