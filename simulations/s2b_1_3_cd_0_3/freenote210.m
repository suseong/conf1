
for k=1:27
    a = mod((k+3),10);
    b = k+3-a;
    fname = ['b2s_1_3__',num2str(b/10),'_',num2str(a)];
    load(fname)
    PPP_{k} = PPP{1};
    RHO_{k} = RHO{1};
    SSS_{k} = SSS{1};
end

%%
xx = []; yy = []; zz = [];

for acc = 1:length(PPP_)
    for num = 1:120
        p = PPP_{acc}(6*(num-1)+1:6*num)/RHO_{acc}(num);
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
            xx{acc}(j,num) = temp(j);      
        end
    end
end

%%
gcf = figure(2);clf;
aa = [1 3 4 6];
bb = [1 1 2 2];
cc = [1.5 1.5 5 5];
for k=1:4
    for j=1:length(PPP_)
        subplot(4,1,k)
        hold on
        plot(linspace(0,6,120),[xx{j}(aa(k),1:120)]);
        axis([0 6 0 cc(k)])
    end
end

%%
% for k=1:27
%     for kk=1:120
%         funnel_s{k}(:,kk) = PPP_{k}(6*(kk-1)+1:6*kk)/RHO_{k}(kk);
%     end
% end

% for k=1:27
%     for kk=1:120
%         funnel_b{k}(:,kk) = PPP_{k}(6*(kk-1)+1:6*kk)/RHO_{k}(kk);
%     end
% end
