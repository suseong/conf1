clear all
close all
% clc

T(1) = 0;
T(2) = 1.0;
T(3) = 1.0;
T(4) = 1.0;
T(5) = 1.0;
T_orig = T;

n = 13;

p1 = [0.0  0.0  1.0];
p2 = [3.0  1.0  1.0];
p3 = [4.5  0.0  1.0];
p4 = [3.0 -1.0  1.0];
p5 = [0.0  0.0  1.0];

p_orig = [p1;p2;p3;p4;p5];
p_ = [p1;p2;p2;p3;p3;p4;p4;p5];

for kkkkk=1:6
    T_ = T_orig*(1+0.2*(kkkkk-1));
    [cost,T_,coeffs_x,coeffs_y,coeffs_z] = calc_min_cost8(T_,p_);
    T_chad = T_(2:end);

%%
for i=1:4
    c_x(:,i) = coeffs_x(n*(i-1)+1:n*i);
    c_y(:,i) = coeffs_y(n*(i-1)+1:n*i);
    c_z(:,i) = coeffs_z(n*(i-1)+1:n*i);    
end

%%

p = p_orig';

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
    
    figure(2)
    subplot(4,3,1)
    hold on
    plot(sum(T(1:k))+T_,pos_x(k,:))
    subplot(4,3,2)
    hold on
    plot(sum(T(1:k))+T_,pos_y(k,:))
    subplot(4,3,3)
    hold on
    plot(sum(T(1:k))+T_,pos_z(k,:))    
    
    subplot(4,3,4)
    hold on
    plot(sum(T(1:k))+T_,vel_x(k,:))
    subplot(4,3,5)
    hold on
    plot(sum(T(1:k))+T_,vel_y(k,:))
    subplot(4,3,6)
    hold on
    plot(sum(T(1:k))+T_,vel_z(k,:))    
    
    subplot(4,3,7)
    hold on
    plot(sum(T(1:k))+T_,acc_x(k,:))
    subplot(4,3,8)
    hold on
    plot(sum(T(1:k))+T_,acc_y(k,:))
    subplot(4,3,9)
    hold on
    plot(sum(T(1:k))+T_,acc_z(k,:))
    
    subplot(4,3,10)
    hold on
    plot(sum(T(1:k))+T_,jerk_x(k,:))
    subplot(4,3,11)
    hold on
    plot(sum(T(1:k))+T_,jerk_y(k,:))
    subplot(4,3,12)
    hold on
    plot(sum(T(1:k))+T_,jerk_z(k,:))    
    
    for kk=1:size(vel_x,2)
       absVel(k,kk) = norm([vel_x(k,kk),vel_y(k,kk),vel_z(k,kk)]);
       absAcc(k,kk) = norm([acc_x(k,kk),acc_y(k,kk),acc_z(k,kk)+9.8]);
    end
    
    figure(3)
    subplot(3,1,1)
    hold on
    plot(sum(T(1:k))+T_,absVel(k,:))
    subplot(3,1,2)
    hold on
    plot(sum(T(1:k))+T_,absAcc(k,:))
    subplot(3,1,3)
    hold on; alpha = 0.35; sinPhi = 0.0523; Delta = 0.5;
    plot(sum(T(1:k))+T_,sinPhi*absAcc(k,:)+alpha*absVel(k,:)+Delta)
    T__ = [T__ T_+sum(T(1:k))];
end

%%
% alpha = 0.35; sinPhi = 0.05; Delta = 0.3;
k=1; disturbance1 = sinPhi*absAcc(k,:)+alpha*absVel(k,:)+Delta;
k=2; disturbance2 = sinPhi*absAcc(k,:)+alpha*absVel(k,:)+Delta;
k=3; disturbance3 = sinPhi*absAcc(k,:)+alpha*absVel(k,:)+Delta;
k=4; disturbance4 = sinPhi*absAcc(k,:)+alpha*absVel(k,:)+Delta;

disturbance = [disturbance1 disturbance2 disturbance3 disturbance4];

ts_d = timeseries(disturbance',T__);
rs_d = resample(ts_d,0:0.03:T__(end));

figure(3)
subplot(3,1,3)
plot([0 rs_d.Time(end)],[0.5 0.5]*1)
plot([0 rs_d.Time(end)],[0.5 0.5]*2)
plot([0 rs_d.Time(end)],[0.5 0.5]*3)
plot([0 rs_d.Time(end)],[0.5 0.5]*4)
plot([0 rs_d.Time(end)],[0.5 0.5]*5)
plot([0 rs_d.Time(end)],[0.5 0.5]*6)

% genCSV(T_chad,c_x,c_y,c_z,['coeffs',num2str(kkkkk),'.csv']);
for k=1:5
    TT(k) = sum(T(1:k));
end

% load funnel_1_3_03.mat
load 'funnel_cd032_deg3.mat'

curRegion = diag(1./([0.05 0.05 0.05 0.1 0.1 0.1]).^2);
funnelIdx = [];
cnt = 1;
for k=1:length(rs_d.Time)
    d = rs_d.Data(k);
%     for kk=1:length(funnel_b)-1
%        if d < 0.1*kk+0.3 && d > 0.1*(kk-1)+0.3
%            d = 0.1*kk+0.3;
%            kk_ = kk;
%        end
%     end

    for kk=1:length(funnel_b)-1
       if d < 0.1*kk+0.4 && d > 0.1*(kk-1)+0.4
           d = 0.1*kk+0.4;
           kk_ = kk;
           break;
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
            
            if  1-1e-3 < norm(chad) && norm(chad) < 1+1e-3
                curRegion = toEllipsoid(funnel_s{kk_}(:,kk+1));
                funnelIdx(1,k) = 1;
                funnelIdx(2,k) = kk_;
                funnelIdx(3,k) = kk+1;
                break;     
            elseif eVal <= 1.00
                curRegion = region;
                funnelIdx(1,k) = 1;
                funnelIdx(2,k) = kk_;
                funnelIdx(3,k) = kk+1;
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

            if 1-1e-3 < norm(chad) && norm(chad) < 1+1e-3
                curRegion = toEllipsoid(funnel_b{kk_}(:,kk+1));
                funnelIdx(1,k) = 2;
                funnelIdx(2,k) = kk_;
                funnelIdx(3,k) = kk+1;
                break;
            elseif eVal <= 1.00
                curRegion = region;
                funnelIdx(1,k) = 2;
                funnelIdx(2,k) = kk_;
                funnelIdx(3,k) = kk+1;
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
    
    if rs_d.Time(k) >= sum(T(1:3))-0.2 && rs_d.Time(k) < sum(T(1:5))+0.2
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
    
       segInfo{kkkkk}{cnt}{1} = [px py pz];
       segInfo{kkkkk}{cnt}{2} = funnelIdx(:,k);
       cnt = cnt+1;
    end
    
    genCSV(T_chad,c_x,c_y,c_z,['coeffs',num2str(kkkkk),'.csv']);

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
    
    figure(101+kkkkk);
    hold on
    surf(fx+px,fy+py,fz+pz-2,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.15,'LineStyle','none');
end

figure(101+kkkkk)
[XX,YY,ZZ] = sphere();
surf(0.2*XX+4,0.2*YY,0.2*ZZ,'FaceColor',[0 0 0],'LineStyle','none');

grid on
axis equal

figure(101+kkkkk)
for k=1:4
    T_ = linspace(0,T(k+1),50);
    p_x(k,:) = calc_trj_(c_x(:,k),T_);
    p_y(k,:) = calc_trj_(c_y(:,k),T_);
    p_z(k,:) = calc_trj_(c_z(:,k),T_); 
    plot3(p_x(k,:),p_y(k,:),p_z(k,:),'-g','linewidth',2);
    axis([-1 6 -3 3 -10 10])
    axis equal
    view(-90,90)
end

% waypoints = {'ok','ok','ob','ob','or'};
% for i=[1 3 5]
%     plot3(p(1,i),p(2,i),p(3,i),waypoints{i},'markersize',8,'linewidth',2)
% end

end




    
    
