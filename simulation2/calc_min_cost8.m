function [cost,T,c_x,c_y,c_z] = calc_min_cost8(T,p)

n = 4;
k = 1/3;
dt = 0.02;
del_t(1,:) = dt*[0 1 -k -k -k];
del_t(2,:) = dt*[0 -k 1 -k -k];
del_t(3,:) = dt*[0 -k -k 1 -k];
del_t(4,:) = dt*[0 -k -k -k 1];

[cost,c_x,c_y,c_z] = calc_cost8(T,p);

cnt = 0;
sw = 0;
convergence = zeros(1,n);

for i=1:500
    cnt = cnt+1;
    if cnt > n
       cnt = 1;
       if (sum(convergence) > 0)
           convergence = zeros(1,n);
       else
           sw = 1;
       end
    end
    
    T_ = T + del_t(cnt,:);
    
    check = min(T_);
    
    if check < 0
        disp('time constraint');
        sw = 1;
    end
    
    if sw == 1
        break;
    end
    [cost_,c_x_,c_y_,c_z_] = calc_cost8(T_,p);
    if(cost_ < cost)
       convergence(cnt) = 1;
       T = T_;
       cost = cost_;
       c_x = c_x_;
       c_y = c_y_;
       c_z = c_z_;
    end
end
end