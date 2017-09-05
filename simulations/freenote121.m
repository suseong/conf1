clear all
close all
clc

for kkkk = 0:8
    
clearvars -except kkkk segInfo
    
T(1) = 0;
T(2) = 0.125;

n = 13;

p1 = [0.0  0.0  1.0];
p2 = [4.5  1.0  2.0];
p3 = [4.5  0.0  2.0];
p4 = [4.5 -1.0  2.0];
p5 = [0.0  0.0  1.0];

p = [p1;p2;p3;p4;p5];
     
v_through = 5 - kkkk*0.3125;
% v_through = 1.0;

v1 = [ 0.0  0.0  0.0];
v2 = [ 0.0 -1.0  0.0]*v_through;
v3 = [ 0.0 -1.0  0.0]*v_through;
v4 = [ 0.0 -1.0  0.0]*v_through;
v5 = [ 0.0  0.0  0.0];

p_ = [p1;p2];
v_ = [v1;v2];

for k=1:20
    T_ = T*(5+k);
    [cost,T__(k,:),~,~,~] = calc_min_cost8(T_,p_,v_);
    costs(k) = cost + 1e6*sum(T__(k,:));
end

figure(1000)
subplot(1,2,1)
plot(costs)

[~,k] = min(costs); 
T_1 =  T__(k,:);
[cost,T_,coeffs_x_1,coeffs_y_1,coeffs_z_1] = calc_min_cost81(T_1,p_,v_);

p_ = [p4;p5];
v_ = [v4;v5];

for k=1:20
    T_ = T*(5+k);
    [cost,T__(k,:),~,~,~] = calc_min_cost8(T_,p_,v_);
    costs(k) = cost + 1e5*sum(T__(k,:));
end

subplot(1,2,2)
plot(costs)

[~,k] = min(costs); 
T_2 =  T__(k,:);
[cost,T_,coeffs_x_2,coeffs_y_2,coeffs_z_2] = calc_min_cost81(T_2,p_,v_);

%%
for i=1:1
    c_x_1(:,i) = coeffs_x_1(n*(i-1)+1:n*i);
    c_y_1(:,i) = coeffs_y_1(n*(i-1)+1:n*i);
    c_z_1(:,i) = coeffs_z_1(n*(i-1)+1:n*i);
    
    c_x_2(:,i) = coeffs_x_2(n*(i-1)+1:n*i);
    c_y_2(:,i) = coeffs_y_2(n*(i-1)+1:n*i);
    c_z_2(:,i) = coeffs_z_2(n*(i-1)+1:n*i);    
end

c_x_3(:,1) = zeros(13,1); c_x_3(1) = p2(1); c_x_3(2) = v2(1);
c_y_3(:,1) = zeros(13,1); c_y_3(1) = p2(2); c_y_3(2) = v2(2); 
c_z_3(:,1) = zeros(13,1); c_z_3(1) = p2(3); c_z_3(2) = v2(3);

c_x_4(:,1) = zeros(13,1); c_x_4(1) = p3(1); c_x_4(2) = v3(1);
c_y_4(:,1) = zeros(13,1); c_y_4(1) = p3(2); c_y_4(2) = v3(2); 
c_z_4(:,1) = zeros(13,1); c_z_4(1) = p3(3); c_z_4(2) = v3(3);

c_x = [c_x_1 c_x_3 c_x_4 c_x_2];
c_y = [c_y_1 c_y_3 c_y_4 c_y_2];
c_z = [c_z_1 c_z_3 c_z_4 c_z_2];

t2 = norm(p3-p2)/norm(v3);
t3 = norm(p4-p3)/norm(v4);

T_ = [T_1 t2 t3 T_2(2:2)];
T_chad = T_(2:end);

%%
p = p';

figure(1)
hold on

T = T_;

for k=1:4
    T_ = linspace(0,T(k+1),50);
    p_x(k,:) = calc_trj_(c_x(:,k),T_);
    p_y(k,:) = calc_trj_(c_y(:,k),T_);
    p_z(k,:) = calc_trj_(c_z(:,k),T_); 
    plot3(p_x(k,:),p_y(k,:),p_z(k,:),'-k');
end

for i=1:5
    plot3(p(1,i),p(2,i),p(3,i),'*','markersize',10)
end

grid on
axis equal

%%
T__ = [];
for k=1:4
    T_ = linspace(0,T(k+1),500);

    pos_x(k,:) = calc_trj_(c_x(:,k),T_);
    pos_y(k,:) = calc_trj_(c_y(:,k),T_);
    pos_z(k,:) = calc_trj_(c_z(:,k),T_); 
    
    vel_x(k,:) = calc_v_(c_x(:,k),T_);
    vel_y(k,:) = calc_v_(c_y(:,k),T_);
    vel_z(k,:) = calc_v_(c_z(:,k),T_); 
    
    acc_x(k,:) = calc_a_(c_x(:,k),T_);
    acc_y(k,:) = calc_a_(c_y(:,k),T_);
    acc_z(k,:) = calc_a_(c_z(:,k),T_);
    
    jerk_x(k,:) = calc_j_(c_x(:,k),T_);
    jerk_y(k,:) = calc_j_(c_y(:,k),T_);
    jerk_z(k,:) = calc_j_(c_z(:,k),T_);
    
%     figure(2)
%     subplot(4,3,1)
%     hold on
%     plot(sum(T(1:k))+T_,pos_x(k,:))
%     subplot(4,3,2)
%     hold on
%     plot(sum(T(1:k))+T_,pos_y(k,:))
%     subplot(4,3,3)
%     hold on
%     plot(sum(T(1:k))+T_,pos_z(k,:))    
%     
%     subplot(4,3,4)
%     hold on
%     plot(sum(T(1:k))+T_,vel_x(k,:))
%     subplot(4,3,5)
%     hold on
%     plot(sum(T(1:k))+T_,vel_y(k,:))
%     subplot(4,3,6)
%     hold on
%     plot(sum(T(1:k))+T_,vel_z(k,:))    
%     
%     subplot(4,3,7)
%     hold on
%     plot(sum(T(1:k))+T_,acc_x(k,:))
%     subplot(4,3,8)
%     hold on
%     plot(sum(T(1:k))+T_,acc_y(k,:))
%     subplot(4,3,9)
%     hold on
%     plot(sum(T(1:k))+T_,acc_z(k,:))
%     
%     subplot(4,3,10)
%     hold on
%     plot(sum(T(1:k))+T_,jerk_x(k,:))
%     subplot(4,3,11)
%     hold on
%     plot(sum(T(1:k))+T_,jerk_y(k,:))
%     subplot(4,3,12)
%     hold on
%     plot(sum(T(1:k))+T_,jerk_z(k,:))    
    
    for kk=1:size(vel_x,2)
       absVel(k,kk) = norm([vel_x(k,kk),vel_y(k,kk),vel_z(k,kk)]);
       absAcc(k,kk) = norm([acc_x(k,kk),acc_y(k,kk),acc_z(k,kk)+9.8]);
    end
    
    figure(3)
    subplot(4,1,1)
    hold on
    plot(sum(T(1:k))+T_,absVel(k,:)+0.01*kkkk)
    subplot(4,1,2)
    hold on
    plot(sum(T(1:k))+T_,absAcc(k,:))
    subplot(4,1,3)
    hold on; alpha = 0.25; sinPhi = 0.03; Delta = 0.2;
    plot(sum(T(1:k))+T_,sinPhi*absAcc(k,:)+alpha*absVel(k,:)+Delta)
    T__ = [T__ T_+sum(T(1:k))];
end

plot([0 T__(end)],[0.5 0.5]*1)
plot([0 T__(end)],[0.5 0.5]*2)
plot([0 T__(end)],[0.5 0.5]*3)
plot([0 T__(end)],[0.5 0.5]*4)
plot([0 T__(end)],[0.5 0.5]*5)

%%
alpha = 0.25; sinPhi = 0.03; Delta = 0.2;
k=1; disturbance1 = sinPhi*absAcc(k,:)+alpha*absVel(k,:)+Delta;
k=2; disturbance2 = sinPhi*absAcc(k,:)+alpha*absVel(k,:)+Delta;
k=3; disturbance3 = sinPhi*absAcc(k,:)+alpha*absVel(k,:)+Delta;
k=4; disturbance4 = sinPhi*absAcc(k,:)+alpha*absVel(k,:)+Delta;
% k=5; disturbance5 = sinPhi*absAcc(k,:)+alpha*absVel(k,:)+Delta;
% k=6; disturbance6 = sinPhi*absAcc(k,:)+alpha*absVel(k,:)+Delta;

% disturbance = [disturbance1 disturbance2 disturbance3 disturbance4 disturbance5 disturbance6];
disturbance = [disturbance1 disturbance2 disturbance3 disturbance4];

ts_d = timeseries(disturbance',T__);
rs_d = resample(ts_d,0:0.03:T__(end));

plot([0 rs_d.Time(end)],[0.5 0.5]*1)
plot([0 rs_d.Time(end)],[0.5 0.5]*2)
plot([0 rs_d.Time(end)],[0.5 0.5]*3)
plot([0 rs_d.Time(end)],[0.5 0.5]*4)
plot([0 rs_d.Time(end)],[0.5 0.5]*5)

%%
for k=1:5
    TT(k) = sum(T(1:k));
end

load funnel_24.mat
curRegion = diag(1./([0.05 0.05 0.05 0.1 0.1 0.1]).^2);
funnelIdx = [];
cnt = 1;
for k=1:length(rs_d.Time)
    d = rs_d.Data(k);
    for kk=1:length(funnel_b)-1
       if d < 0.1*kk+0.3 && d > 0.1*(kk-1)+0.3
           d = 0.1*kk+0.2;
           kk_ = kk;
       end
    end
    u_bound = toEllipsoid(funnel_s{kk_}(:,120));
    % check whether curRegion is smaller than u_bound or not
    L = chol(curRegion)';
    invL = inv(L);
    eVal = max(eig(invL*u_bound*invL'));
    
    if eVal < 1
        smallerThanUBound = 1;
    else
        smallerThanUBound = 0;
    end
    
    if smallerThanUBound
        for kk=10:120
            region = toEllipsoid(funnel_s{kk_}(:,kk));
            L = chol(curRegion)';
            invL = inv(L);
            chad = invL*region*invL';
            eVal = max(eig(chad));
            
            if  1-1e-2 < norm(chad) && norm(chad) < 1+1e-2
                curRegion = toEllipsoid(funnel_s{kk_}(:,kk+1));
                funnelIdx(1,k) = 1;
                funnelIdx(2,k) = kk_;
                funnelIdx(3,k) = kk+1;
                break;     
            elseif eVal <= 1.01
                curRegion = region;
                funnelIdx(1,k) = 1;
                funnelIdx(2,k) = kk_;
                funnelIdx(3,k) = kk;
                break;
            elseif eVal > 1
            else
                keyboard;
            end
        end
    else
        for kkk=1:110
            kk = 121-kkk;
            if kk >= 118
                kk = 118;
            end
            region = toEllipsoid(funnel_b{kk_}(:,kk));
            L = chol(curRegion)';
            invL = inv(L);
            chad = invL*region*invL';
            eVal = max(eig(chad));

            if 1-1e-2 < norm(chad) && norm(chad) < 1+1e-2
                curRegion = toEllipsoid(funnel_b{kk_}(:,kk+1));
                funnelIdx(1,k) = 2;
                funnelIdx(2,k) = kk_;
                funnelIdx(3,k) = kk+1;
                break;
            elseif eVal <= 1.01
                curRegion = region;
                funnelIdx(1,k) = 2;
                funnelIdx(2,k) = kk_;
                funnelIdx(3,k) = kk;
                break;
            elseif eVal > 1
            else
                keyboard;
            end            
        end
    end
    if length(funnelIdx) < k
        keyboard
    end
    
    if rs_d.Time(k) >= sum(T(1:2))-0.2 && rs_d.Time(k) < sum(T(1:4))+0.2
        t = 0.03*k;
        
        for kk=1:4
            if t > TT(kk) && t <= TT(kk+1)
                idx = kk;
            end
        end
        
        t = t-TT(idx);
                
        px = calc_trj_(c_x(:,idx),t);
        py = calc_trj_(c_y(:,idx),t);
        pz = calc_trj_(c_z(:,idx),t);
    
       segInfo{kkkk+1}{cnt}{1} = [px py pz];
       segInfo{kkkk+1}{cnt}{2} = funnelIdx(:,k);
       cnt = cnt+1;
    end
    
    genCSV(T_chad,c_x,c_y,c_z,['coeffs',num2str(kkkk),'.csv']);

end
 
%%
for k=1:5
    TT(k) = sum(T(1:k));
end

% for k=1:length(funnelIdx)
for kkk=1:length(funnelIdx)
    k = kkk;
    if funnelIdx(1,k) == 1
        chad = toEllipsoid(funnel_s{funnelIdx(2,k)}(:,funnelIdx(3,k)));
        Pp = chad(1:3,1:3); Ppv = chad(1:3,4:6); Pv = chad(4:6,4:6);
        Px = Pp-Ppv^2*pinv(Pv);
        [fx,fy,fz] = sphere(10);
        fx = fx/sqrt(Px(1,1));
        fy = fy/sqrt(Px(2,2));
        fz = fz/sqrt(Px(3,3));
    else
        chad = toEllipsoid(funnel_b{funnelIdx(2,k)}(:,funnelIdx(3,k)));
        Pp = chad(1:3,1:3); Ppv = chad(1:3,4:6); Pv = chad(4:6,4:6);
        Px = Pp-Ppv^2*pinv(Pv);
        [fx,fy,fz] = sphere(10);
        fx = fx/sqrt(Px(1,1));
        fy = fy/sqrt(Px(2,2));
        fz = fz/sqrt(Px(3,3));
    end
    
    figure(3);
%     subplot(4,1,4)
%     hold on
%     surf(fx+0.3*k,fy,fz)
    
    t = 0.03*k;
    
    for kk=1:4
       if t > TT(kk) && t <= TT(kk+1)
           idx = kk;
       end
    end
    
    t = t-TT(idx);
    
    px = calc_trj_(c_x(:,idx),t);
    py = calc_trj_(c_y(:,idx),t);
    pz = calc_trj_(c_z(:,idx),t); 
    
    figure(101+kkkk);
    hold on
    surf(fx+px,fy+py,fz+pz,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.1,'LineStyle','none');
end

grid on
axis equal

figure(101+kkkk)
for k=1:4
    T_ = linspace(0,T(k+1),50);
    p_x(k,:) = calc_trj_(c_x(:,k),T_);
    p_y(k,:) = calc_trj_(c_y(:,k),T_);
    p_z(k,:) = calc_trj_(c_z(:,k),T_); 
    plot3(p_x(k,:),p_y(k,:),p_z(k,:),'-k');
end

for i=1:5
    plot3(p(1,i),p(2,i),p(3,i),'*','markersize',10)
end

figure(3);
subplot(4,1,4)
grid on
axis equal

end

% genCSV(T_chad,c_x,c_y,c_z,['coeffs',num2str(kkkk),'.csv']);

%%
for kkkk=1:length(segInfo)
    for kkk=1:length(segInfo{kkkk})
        funnelIdx = segInfo{kkkk}{kkk}{2};
        pos{kkkk}(kkk,:) = segInfo{kkkk}{kkk}{1};
        if funnelIdx(1) == 1
            chad = toEllipsoid(funnel_s{funnelIdx(2)}(:,funnelIdx(3)));
        else
            chad = toEllipsoid(funnel_b{funnelIdx(2)}(:,funnelIdx(3)));
        end
        Pp = chad(1:3,1:3); Ppv = chad(1:3,4:6); Pv = chad(4:6,4:6);
        Px = Pp-Ppv^2*pinv(Pv);
        funnelSize{kkkk}(:,kkk) = 1./sqrt(diag(Px));
    end
end

%%
figure(10001)
for k=1:length(segInfo)
    subplot(2,1,1)
    hold on
    plot(-pos{k}(:,2),funnelSize{k}(1,:),'-')
    subplot(2,1,2)
    hold on
    plot(-pos{k}(:,2),funnelSize{k}(3,:),'-')    
end

subplot(2,1,1)
grid on
axis([-1 1 0.0 0.6])
subplot(2,1,2)
grid on
axis([-1 1 0.0 0.5])

