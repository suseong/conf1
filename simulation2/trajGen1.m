clear all
close all

syms t real
syms wj ws real

T = [];
A = []; %1
B = []; %2
C = []; %3
D = []; %4

n = 13;

for i=0:n-1
    T = [T t^i];
    A = [A sym(['a',num2str(i)],'real')];
    B = [B sym(['b',num2str(i)],'real')];
    C = [C sym(['c',num2str(i)],'real')];
    D = [D sym(['d',num2str(i)],'real')];
end

trj_a = A*T';
trj_b = B*T';
trj_c = C*T';
trj_d = D*T';

d1trj_a = diff(trj_a,t);
d2trj_a = diff(d1trj_a,t);
d3trj_a = diff(d2trj_a,t);
d4trj_a = diff(d3trj_a,t);

d1trj_b = diff(trj_b,t);
d2trj_b = diff(d1trj_b,t);
d3trj_b = diff(d2trj_b,t);
d4trj_b = diff(d3trj_b,t);

d1trj_c = diff(trj_c,t);
d2trj_c = diff(d1trj_c,t);
d3trj_c = diff(d2trj_c,t);
d4trj_c = diff(d3trj_c,t);

d1trj_d = diff(trj_d,t);
d2trj_d = diff(d1trj_d,t);
d3trj_d = diff(d2trj_d,t);
d4trj_d = diff(d3trj_d,t);

sq_d4trj_a = d4trj_a*d4trj_a;
sq_d4trj_b = d4trj_b*d4trj_b;
sq_d4trj_c = d4trj_c*d4trj_c;
sq_d4trj_d = d4trj_d*d4trj_d;

chad_a = int(sq_d4trj_a);
chad_b = int(sq_d4trj_b);
chad_c = int(sq_d4trj_c);
chad_d = int(sq_d4trj_d);
 
syms t0 t1 t2 t3 t4 real;

int_sq_d3trj(1) = subs(chad_a,t,t1) - subs(chad_a,t,0);
int_sq_d3trj(2) = subs(chad_b,t,t2) - subs(chad_b,t,0);
int_sq_d3trj(3) = subs(chad_c,t,t3) - subs(chad_c,t,0);
int_sq_d3trj(4) = subs(chad_d,t,t4) - subs(chad_d,t,0);

HH = sym(zeros(4*n));

a_h = ['a','b','c','d'];

i=1;
H_ = sym(zeros(n,n));
for j=1:n
    for k=1:n
        chad = sym([a_h(i),num2str(j-1)]);
        aa = diff(int_sq_d3trj(i),chad);
        chad = sym([a_h(i),num2str(k-1)]);
        bb = diff(aa,chad);
        H_(j,k) = H_(j,k) + 0.5*bb;
    end
end

% matlabFunction(H_,'file','calc_H12','Vars',{wj,ws,[t0 t1]})
matlabFunction(H_,'file','calc_H12','Vars',{[t0 t1]})

%%
const = sym([]);
const(1) = subs(trj_a,t,0);  % p(1);
const(2) = subs(trj_a,t,t1);  % p(2);
const(3) = subs(trj_b,t,0);  % p(2);
const(4) = subs(trj_b,t,t2);  % p(3);
const(5) = subs(trj_c,t,0);  % p(3);
const(6) = subs(trj_c,t,t3);  % p(4);
const(7) = subs(trj_d,t,0);  % p(4);
const(8) = subs(trj_d,t,t4);  % p(5);

% const(9)  = subs(d1trj_a,t,0); % v(1);
% const(10) = subs(d1trj_a,t,t1); % v(2);
% const(11) = subs(d1trj_b,t,0); % v(2);
% const(12) = subs(d1trj_b,t,t2); % v(3);
% const(13) = subs(d1trj_c,t,0); % v(3);
% const(14) = subs(d1trj_c,t,t3); % v(4);
% const(15) = subs(d1trj_d,t,0); % v(4);
% const(16) = subs(d1trj_d,t,t4); % v(5);

const(9)  = subs(d1trj_a,t,0);
const(10) = subs(d1trj_a,t,t1) - subs(d1trj_b,t,0);
const(11) = subs(d1trj_b,t,t2) - subs(d1trj_c,t,0);
const(12) = subs(d1trj_c,t,t3) - subs(d1trj_d,t,0);
const(13) = subs(d1trj_d,t,t4);

const(14) = subs(d2trj_a,t,0);
const(15) = subs(d2trj_a,t,t1) - subs(d2trj_b,t,0);
const(16) = subs(d2trj_b,t,t2) - subs(d2trj_c,t,0);
const(17) = subs(d2trj_c,t,t3) - subs(d2trj_d,t,0);
const(18) = subs(d2trj_d,t,t4);

const(19) = subs(d3trj_a,t,t0);
const(20) = subs(d3trj_a,t,t1) - subs(d3trj_b,t,0);
const(21) = subs(d3trj_b,t,t2) - subs(d3trj_c,t,0);
const(22) = subs(d3trj_c,t,t3) - subs(d3trj_d,t,0);
const(23) = subs(d3trj_d,t,t4);

const(24) = subs(d4trj_a,t,t0);
const(25) = subs(d4trj_a,t,t1) - subs(d4trj_b,t,0);
const(26) = subs(d4trj_b,t,t2) - subs(d4trj_c,t,0);
const(27) = subs(d4trj_c,t,t3) - subs(d4trj_d,t,0);
const(28) = subs(d4trj_d,t,t4);

% const(18) = subs(d2trj_a - d2trj_b,t,t1);
% const(19) = subs(d2trj_b - d2trj_c,t,t2);
% const(20) = subs(d2trj_c - d2trj_d,t,t3);
% const(21) = subs(d2trj_d,t,t4);
% 
% const(22) = subs(d3trj_a,t,t0);
% const(23) = subs(d3trj_a - d3trj_b,t,t1);
% const(24) = subs(d3trj_b - d3trj_c,t,t2);
% const(25) = subs(d3trj_c - d3trj_d,t,t3);
% const(26) = subs(d3trj_d,t,t4);
% 
% const(27) = subs(d4trj_a,t,t0);
% const(28) = subs(d4trj_a - d4trj_b,t,t1);
% const(29) = subs(d4trj_b - d4trj_c,t,t2);
% const(30) = subs(d4trj_c - d4trj_d,t,t3);
% const(31) = subs(d4trj_d,t,t4);
 
variables = [A B C D];

AA = sym(zeros(size(const,2),size(variables,2)));

for i=1:size(const,2)
    for j=1:size(variables,2)
        AA(i,j) = diff(const(i),variables(j));
    end
end

matlabFunction(AA,'file','calc_AA12_','Vars',{[t0 t1 t2 t3 t4]})

%%
% p(1,:) = [2   1 1];
% p(2,:) = [1   0 0];
% p(3,:) = [0.3 0 0];
% p(4,:) = [0.2 0 0];
% p(5,:) = [0.1 0 0];
% p(6,:) = [0   0 0];
% p(7,:) = [0   0 0];
% p(8,:) = [0   0 0];
% p(9,:) = [0   0 0];
% 
% %% traj constraint
% const = sym([]);
% const(1) = subs(trj_a,t,t0);  % p(1);
% const(2) = subs(trj_a,t,t1);  % p(2);
% const(3) = subs(trj_b,t,t1);  % p(2);
% const(4) = subs(trj_b,t,t2);  % p(3);
% const(5) = subs(trj_c,t,t2);  % p(3);
% const(6) = subs(trj_c,t,t3);  % p(4);
% const(7) = subs(trj_d,t,t3);  % p(4);
% const(8) = subs(trj_d,t,t4);  % p(5);
% const(9) = subs(trj_e,t,t4);  % p(5);
% const(10) = subs(trj_e,t,t5); % p(6);
% const(11) = subs(trj_f,t,t5); % p(6);
% const(12) = subs(trj_f,t,t6); % p(7);
% const(13) = subs(trj_g,t,t6); % p(7);
% const(14) = subs(trj_g,t,t7); % p(8);
% const(15) = subs(trj_h,t,t7); % p(8);
% const(16) = subs(trj_h,t,t8); % p(9);
% 
% const(17) = subs(d1trj_a,t,t0); % v(1);
% const(18) = subs(d1trj_a,t,t1); % v(2);
% const(19) = subs(d1trj_b,t,t1); % v(2);
% const(20) = subs(d1trj_b,t,t2); % v(3);
% const(21) = subs(d1trj_c,t,t2); % v(3);
% const(22) = subs(d1trj_c,t,t3); % v(4);
% const(23) = subs(d1trj_d,t,t3); % v(4);
% const(24) = subs(d1trj_d,t,t4); % v(5);
% const(25) = subs(d1trj_e,t,t4); % v(5);
% const(26) = subs(d1trj_e,t,t5); % v(6);
% const(27) = subs(d1trj_f,t,t5); % v(6);
% const(28) = subs(d1trj_f,t,t6); % v(7);
% const(29) = subs(d1trj_g,t,t6); % v(7);
% const(30) = subs(d1trj_g,t,t7); % v(8);
% const(31) = subs(d1trj_h,t,t7); % v(8);
% const(32) = subs(d1trj_h,t,t8); % v(9);
% 
% const(33) = subs(d2trj_a,t,t0);
% const(34) = subs(d2trj_a - d2trj_b,t,t1);
% const(35) = subs(d2trj_b - d2trj_c,t,t2);
% const(36) = subs(d2trj_c - d2trj_d,t,t3);
% const(37) = subs(d2trj_d - d2trj_e,t,t4);
% const(38) = subs(d2trj_e - d2trj_f,t,t5);
% const(39) = subs(d2trj_f - d2trj_g,t,t6);
% const(40) = subs(d2trj_g - d2trj_h,t,t7);
% const(41) = subs(d2trj_h,t,t8);
% 
% const(42) = subs(d3trj_a,t,t0);
% const(43) = subs(d3trj_a - d3trj_b,t,t1);
% const(44) = subs(d3trj_b - d3trj_c,t,t2);
% const(45) = subs(d3trj_c - d3trj_d,t,t3);
% const(46) = subs(d3trj_d - d3trj_e,t,t4);
% const(47) = subs(d3trj_e - d3trj_f,t,t5);
% const(48) = subs(d3trj_f - d3trj_g,t,t6);
% const(49) = subs(d3trj_g - d3trj_h,t,t7);
% const(50) = subs(d3trj_h,t,t8);
% 
% variables = [A B C D E F G H];
% 
% AA = sym(zeros(size(const,2),size(variables,2)));
% 
% for i=1:size(const,2)
%     for j=1:size(variables,2)
%         AA(i,j) = diff(const(i),variables(j));
%     end
% end

%%

% %%
% l = 1;
% aa = [p(1,l);p(2,l);p(2,l);p(3,l);p(3,l);p(4,l);p(4,l);p(5,l);p(5,l);p(6,l);zeros(size(const,2)-10,1)];
% 
% x = quadprog(H,zeros(1,size(H,1)),zeros(size(H)),zeros(size(H,1),1),double(AA),aa);
% 
% %%
% i=1; a_c = x(n*(i-1)+1:n*i);
% i=2; b_c = x(n*(i-1)+1:n*i);
% i=3; c_c = x(n*(i-1)+1:n*i);
% i=4; d_c = x(n*(i-1)+1:n*i);
% i=5; e_c = x(n*(i-1)+1:n*i);
% 
% %%
% figure(1)
% hold on
% Ta = t0:0.02:t1;
% Tb = t1:0.02:t2;
% Tc = t2:0.02:t3;
% Td = t3:0.02:t4;
% Te = t4:0.02:t5;
% 
% p1 = calc_trj(a_c,Ta);
% p2 = calc_trj(b_c,Tb);
% p3 = calc_trj(c_c,Tc);
% p4 = calc_trj(d_c,Td);
% p5 = calc_trj(e_c,Te);
% 
% plot(Ta,p1);
% plot(Tb,p2);
% plot(Tc,p3);
% plot(Td,p4);
% plot(Te,p5);
% 
% 
% %%
% double(subs(trj_a,[A t],[a_c' t0]))
% double(subs(trj_a,[A t],[a_c' t1]))
% double(subs(trj_b,[B t],[b_c' t1]))
% double(subs(trj_b,[B t],[b_c' t2]))
% double(subs(trj_c,[C t],[c_c' t2]))
% double(subs(trj_c,[C t],[c_c' t3]))
% double(subs(trj_d,[D t],[d_c' t3]))
% double(subs(trj_d,[D t],[d_c' t4]))
% double(subs(trj_e,[E t],[e_c' t4]))
% double(subs(trj_e,[E t],[e_c' t5]))
% 
% %%
% double(subs(d1trj_a,[A t],[a_c' t0]))
% double(subs(d1trj_a,[A t],[a_c' t1]))
% double(subs(d1trj_b,[B t],[b_c' t1]))
% double(subs(d1trj_b,[B t],[b_c' t2]))
% double(subs(d1trj_c,[C t],[c_c' t2]))
% double(subs(d1trj_c,[C t],[c_c' t3]))
% double(subs(d1trj_d,[D t],[d_c' t3]))
% double(subs(d1trj_d,[D t],[d_c' t4]))
% double(subs(d1trj_e,[E t],[e_c' t4]))
% double(subs(d1trj_e,[E t],[e_c' t5]))
% 
% %%
% double(subs(d2trj_a,[A t],[a_c' t0]))
% double(subs(d2trj_a,[A t],[a_c' t1]))
% double(subs(d2trj_b,[B t],[b_c' t1]))
% double(subs(d2trj_b,[B t],[b_c' t2]))
% double(subs(d2trj_c,[C t],[c_c' t2]))
% double(subs(d2trj_c,[C t],[c_c' t3]))
% double(subs(d2trj_d,[D t],[d_c' t3]))
% double(subs(d2trj_d,[D t],[d_c' t4]))
% double(subs(d2trj_e,[E t],[e_c' t4]))
% double(subs(d2trj_e,[E t],[e_c' t5]))
% 
% %%
% double(subs(d3trj_a,[A t],[a_c' t0]))
% double(subs(d3trj_a,[A t],[a_c' t1]))
% double(subs(d3trj_b,[B t],[b_c' t1]))
% double(subs(d3trj_b,[B t],[b_c' t2]))
% double(subs(d3trj_c,[C t],[c_c' t2]))
% double(subs(d3trj_c,[C t],[c_c' t3]))
% double(subs(d3trj_d,[D t],[d_c' t3]))
% double(subs(d3trj_d,[D t],[d_c' t4]))
% double(subs(d3trj_e,[E t],[e_c' t4]))
% double(subs(d3trj_e,[E t],[e_c' t5]))
% 
% %%
% 
% 
% 
% 
