load funnel_cd032_deg3

xx = [];

for acc_ = 1:6
    acc = acc_*5-3;
    for num = 1:120
        p = funnel_s{acc}(:,num);
        P = [p(1)   0    0  p(3)   0    0;
            0  p(1)   0    0  p(3)   0;
            0    0  p(2)   0    0  p(4);
            p(3)   0    0  p(5)   0    0;
            0  p(3)   0    0  p(5)   0;
            0    0  p(4)   0    0  p(6)];
        
        invP = inv(P);
        for j = 1:6
            a = zeros(6,1); a(j) = 1;
            lambda = sqrt(a'*invP*a);
            temp = 1/lambda*invP*a;
            xx{acc_}(j,num) = temp(j);      
        end        
    end
end

%%
gcf = figure(1);clf;
aa = [1 3 4 6];
bb = [1 1 2 2];
cc = [1.5 1.5 5 5];
lineChad = {'-k','--k','-b','--b','-g','--g','r'};
for k=1:4
for j=1:6
    subplot(4,1,k)
    hold on
    plot(linspace(0,6,120)-0.05,[xx{j}(aa(k),1:120)],lineChad{j});
%     plot(linspace(0,6,120)-0.05,[xx{j}(aa(k),1:120)]);
    axis([0 5 0 cc(k)])
end
end

%%
gcf = figure(1);

subplot(4,1,1)
yticks([0 0.5 1 1.5])
xticks([0 1 2 3 4 5 6])
ylabel('$e_{p,xy}$ (m)','Interpreter','latex','FontSize',12)
grid on

subplot(4,1,2)
yticks([0 0.5 1 1.5])
xticks([0 1 2 3 4 5 6])
ylabel('$e_{p,z}$ (m)','Interpreter','latex','FontSize',12)
grid on

subplot(4,1,3)
yticks([0 2.5 5])
xticks([0 1 2 3 4 5 6])
ylabel('$e_{v,xy}$ (m/s)','Interpreter','latex','FontSize',12)
grid on

subplot(4,1,4)
yticks([0 2.5 5])
xticks([0 1 2 3 4 5 6]) 
ylabel('$e_{v,z}$ (m/s)','Interpreter','latex','FontSize',12)
grid on

xlabel('time (sec)','Interpreter','latex','FontSize',12)
 
set(gcf, 'Position', [100, 100, 450, 400]);

%%
xx = [];

for acc_ = 1:6
    acc = acc_*5-3;
    for num = 1:120
        p = funnel_b{acc}(:,num);
        P = [p(1)   0    0  p(3)   0    0;
            0  p(1)   0    0  p(3)   0;
            0    0  p(2)   0    0  p(4);
            p(3)   0    0  p(5)   0    0;
            0  p(3)   0    0  p(5)   0;
            0    0  p(4)   0    0  p(6)];
        
        invP = inv(P);
        for j = 1:6
            a = zeros(6,1); a(j) = 1;
            lambda = sqrt(a'*invP*a);
            temp = 1/lambda*invP*a;
            xx{acc_}(j,num) = temp(j);      
        end        
    end
end

%%
gcf = figure(2);clf;
aa = [1 3 4 6];
bb = [1 1 2 2];
cc = [1.5 1.5 5 5];
lineChad = {'-k','--k','-b','--b','-g','--g','r'};
for k=1:4
for j=1:6
    subplot(4,1,k)
    hold on
    plot(linspace(0,6,120)-0.2,[xx{j}(aa(k),1:120)],lineChad{j});
%     plot(linspace(0,6,120)-0.20,[xx{j}(aa(k),1:120)]);
    axis([0 5 0 cc(k)])
end
end

%%
gcf = figure(2);

subplot(4,1,1)
yticks([0 0.5 1 1.5])
xticks([0 1 2 3 4 5 6])
ylabel('$e_{p,xy}$ (m)','Interpreter','latex','FontSize',12)
grid on

subplot(4,1,2)
yticks([0 0.5 1 1.5])
xticks([0 1 2 3 4 5 6])
ylabel('$e_{p,z}$ (m)','Interpreter','latex','FontSize',12)
grid on

subplot(4,1,3)
yticks([0 2.5 5])
xticks([0 1 2 3 4 5 6])
ylabel('$e_{v,xy}$ (m/s)','Interpreter','latex','FontSize',12)
grid on

subplot(4,1,4)
yticks([0 2.5 5])
xticks([0 1 2 3 4 5 6]) 
ylabel('$e_{v,z}$ (m/s)','Interpreter','latex','FontSize',12)
grid on

xlabel('time (sec)','Interpreter','latex','FontSize',12)
 
set(gcf, 'Position', [100, 100, 450, 400]);
