function [cost,c_x,c_y,c_z] = calc_cost8(T,p,v)

n = 2;
order = 13;
H = zeros(order*n);

for i=1:n
   H(order*(i-1)+1:order*i,order*(i-1)+1:order*i) = calc_H12([T(i) T(i+1)]); 
end

AA = calc_AA12_(T);
nullAA = null(AA);

AA_ = [AA;nullAA'];

invAA_ = inv(AA_);
H_ = invAA_'*H*invAA_;

a_length = size(p,1) + size(v,1) + 3*(n+1);
b_length = n*order - a_length;

H_a = H_(1:a_length,1:a_length);
H_b = H_(a_length+1:end,a_length+1:end);
H_ab = H_(1:a_length,a_length+1:end);

for k = 1:3
    a = [p(:,k);v(:,k);zeros(3*(n+1),1)];
    b_ = -inv(H_b)*H_ab'*a;
    c_(:,k) = invAA_*[a;b_];
%     cost_(k) = [a;b_]'*invAA_'*H*invAA_*[a;b_];
    cost_(k) = c_(:,k)'*H*c_(:,k);
    
end

c_x = c_(:,1);
c_y = c_(:,2);
c_z = c_(:,3);

cost = sum(cost_);

end