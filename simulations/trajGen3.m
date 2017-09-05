clear all
close all

syms t real

T = [];
A = []; %1
B = []; %2
C = []; %3
D = []; %4
E = []; %5
F = []; %6
G = []; %7

n = 13;

for i=0:n-1
    T = [T t^i];
    A = [A sym(['a',num2str(i)],'real')];
    B = [B sym(['b',num2str(i)],'real')];
    C = [C sym(['c',num2str(i)],'real')];
    D = [D sym(['d',num2str(i)],'real')];
    E = [E sym(['e',num2str(i)],'real')];
    F = [F sym(['f',num2str(i)],'real')];
%     G = [G sym(['g',num2str(i)],'real')];
end

trj_a = A*T';
trj_b = B*T';
trj_c = C*T';
trj_d = D*T';
trj_e = E*T';
trj_f = F*T';
% trj_g = G*T';

d1trj_a = diff(  trj_a,t);
d2trj_a = diff(d1trj_a,t);
d3trj_a = diff(d2trj_a,t);
d4trj_a = diff(d3trj_a,t);

d1trj_b = diff(  trj_b,t);
d2trj_b = diff(d1trj_b,t);
d3trj_b = diff(d2trj_b,t);
d4trj_b = diff(d3trj_b,t);

d1trj_c = diff(  trj_c,t);
d2trj_c = diff(d1trj_c,t);
d3trj_c = diff(d2trj_c,t);
d4trj_c = diff(d3trj_c,t);

d1trj_d = diff(  trj_d,t);
d2trj_d = diff(d1trj_d,t);
d3trj_d = diff(d2trj_d,t);
d4trj_d = diff(d3trj_d,t);

d1trj_e = diff(  trj_e,t);
d2trj_e = diff(d1trj_e,t);
d3trj_e = diff(d2trj_e,t);
d4trj_e = diff(d3trj_e,t);

d1trj_f = diff(  trj_f,t);
d2trj_f = diff(d1trj_f,t);
d3trj_f = diff(d2trj_f,t);
d4trj_f = diff(d3trj_f,t);

% d1trj_g = diff(  trj_g,t);
% d2trj_g = diff(d1trj_g,t);
% d3trj_g = diff(d2trj_g,t);
% d4trj_g = diff(d3trj_g,t);

% sq_d4trj_a = d4trj_a*d4trj_a;
% sq_d4trj_b = d4trj_b*d4trj_b;
% sq_d4trj_c = d4trj_c*d4trj_c;
% sq_d4trj_d = d4trj_d*d4trj_d;
% sq_d4trj_e = d4trj_e*d4trj_e;
% sq_d4trj_f = d4trj_f*d4trj_f;
% sq_d4trj_g = d4trj_g*d4trj_g;
% 
% chad_a = int(sq_d4trj_a);
% chad_b = int(sq_d4trj_b);
% chad_c = int(sq_d4trj_c);
% chad_d = int(sq_d4trj_d);
% chad_e = int(sq_d4trj_e);
% chad_f = int(sq_d4trj_f);
% chad_g = int(sq_d4trj_g);

sq_d4trj_a = d4trj_a*d4trj_a;
sq_d4trj_b = d4trj_b*d4trj_b;
sq_d4trj_c = d4trj_c*d4trj_c;
sq_d4trj_d = d4trj_d*d4trj_d;
sq_d4trj_e = d4trj_e*d4trj_e;
sq_d4trj_f = d4trj_f*d4trj_f;
% sq_d3trj_g = d3trj_g*d3trj_g;

chad_a = int(sq_d4trj_a);
chad_b = int(sq_d4trj_b);
chad_c = int(sq_d4trj_c);
chad_d = int(sq_d4trj_d);
chad_e = int(sq_d4trj_e);
chad_f = int(sq_d4trj_f);
% chad_g = int(sq_d4trj_g);

syms t0 t1 t2 t3 t4 t5 t6 t7 t8 real;

int_sq_d3trj(1) = subs(chad_a,t,t1) - subs(chad_a,t,0);
int_sq_d3trj(2) = subs(chad_b,t,t2) - subs(chad_b,t,0);
int_sq_d3trj(3) = subs(chad_c,t,t3) - subs(chad_c,t,0);
int_sq_d3trj(4) = subs(chad_d,t,t4) - subs(chad_d,t,0);
int_sq_d3trj(5) = subs(chad_e,t,t5) - subs(chad_e,t,0);
int_sq_d3trj(6) = subs(chad_f,t,t6) - subs(chad_f,t,0);
% int_sq_d3trj(7) = subs(chad_g,t,t7) - subs(chad_g,t,t6);

a_h = ['a','b','c','d','e','f','g'];

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

matlabFunction(H_,'file','calc_H12','Vars',{[t0 t1]})

%% traj constraint
const = sym([]);
const(1)  = subs(trj_a,t,0);  % p(1);
const(2)  = subs(trj_a,t,t1);  % p(2);
const(3)  = subs(trj_b,t,0);  % p(2);
const(4)  = subs(trj_b,t,t2);  % p(3);
const(5)  = subs(trj_c,t,0);  % p(3);
const(6)  = subs(trj_c,t,t3);  % p(4);
const(7)  = subs(trj_d,t,0);  % p(4);
const(8)  = subs(trj_d,t,t4);  % p(5);
const(9)  = subs(trj_e,t,0);  % p(5);
const(10) = subs(trj_e,t,t5); % p(6);
const(11) = subs(trj_f,t,0); % p(6);
const(12) = subs(trj_f,t,t6); % p(7);
% const(13) = subs(trj_g,t,t6); % p(7);
% const(14) = subs(trj_g,t,t7); % p(8);

const(13) = subs(d1trj_a,t,0); % v(1);
const(14) = subs(d1trj_a,t,t1); % v(2);
const(15) = subs(d1trj_b,t,0); % v(2);
const(16) = subs(d1trj_b,t,t2); % v(3);
const(17) = subs(d1trj_c,t,0); % v(3);
const(18) = subs(d1trj_c,t,t3); % v(4);
const(19) = subs(d1trj_d,t,0); % v(4);
const(20) = subs(d1trj_d,t,t4); % v(5);
const(21) = subs(d1trj_e,t,0); % v(5);
const(22) = subs(d1trj_e,t,t5); % v(6);
const(23) = subs(d1trj_f,t,0); % v(6);
const(24) = subs(d1trj_f,t,t6); % v(7);
% const(27) = subs(d1trj_g,t,t6); % v(7);
% const(28) = subs(d1trj_g,t,t7); % v(8);

const(25) = subs(d2trj_a,t,0);
const(26) = subs(d2trj_a,t,t1) - subs(d2trj_b,t,0);
const(27) = subs(d2trj_b,t,t2) - subs(d2trj_c,t,0);
const(28) = subs(d2trj_c,t,t3) - subs(d2trj_d,t,0);
const(29) = subs(d2trj_d,t,t4) - subs(d2trj_e,t,0);
const(30) = subs(d2trj_e,t,t5) - subs(d2trj_f,t,0);
const(31) = subs(d2trj_f,t,t6);
% const(32) = subs(d2trj_g,t,t7);

const(32) = subs(d3trj_a,t,0);
const(33) = subs(d3trj_a,t,t1) - subs(d3trj_b,t,0);
const(34) = subs(d3trj_b,t,t2) - subs(d3trj_c,t,0);
const(35) = subs(d3trj_c,t,t3) - subs(d3trj_d,t,0);
const(36) = subs(d3trj_d,t,t4) - subs(d3trj_e,t,0);
const(37) = subs(d3trj_e,t,t5) - subs(d3trj_f,t,0);
const(38) = subs(d3trj_f,t,t6);
% const(44) = subs(d3trj_g,t,t7);

const(39) = subs(d4trj_a,t,0);
const(40) = subs(d4trj_a,t,t1) - subs(d4trj_b,t,0);
const(41) = subs(d4trj_b,t,t2) - subs(d4trj_c,t,0);
const(42) = subs(d4trj_c,t,t3) - subs(d4trj_d,t,0);
const(43) = subs(d4trj_d,t,t4) - subs(d4trj_e,t,0);
const(44) = subs(d4trj_e,t,t5) - subs(d4trj_f,t,0);
const(45) = subs(d4trj_f,t,t6);
% const(52) = subs(d4trj_g,t,t7);

variables = [A B C D E F];

AA = sym(zeros(size(const,2),size(variables,2)));

for i=1:size(const,2)
    for j=1:size(variables,2)
        AA(i,j) = diff(const(i),variables(j));
    end
end

matlabFunction(AA,'file','calc_AA12_','Vars',{[t0 t1 t2 t3 t4 t5 t6]})
