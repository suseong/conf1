clear all
close all

syms t real

T = [];
A = []; %1
B = []; %2

n = 13;

for i=0:n-1
    T = [T t^i];
    A = [A sym(['a',num2str(i)],'real')];
end

trj_a = A*T';

d1trj_a = diff(  trj_a,t);
d2trj_a = diff(d1trj_a,t);
d3trj_a = diff(d2trj_a,t);
d4trj_a = diff(d3trj_a,t);

sq_d4trj_a = d4trj_a*d4trj_a;

chad_a = int(sq_d4trj_a);

syms t0 t1 real;

int_sq_d4trj(1) = subs(chad_a,t,t1) - subs(chad_a,t,0);

a_h = ['a'];

i=1;
H_ = sym(zeros(n,n));
for j=1:n
    for k=1:n
        chad = sym([a_h(i),num2str(j-1)]);
        aa = diff(int_sq_d4trj(i),chad);
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

const(3) = subs(d1trj_a,t,0); % v(1);
const(4) = subs(d1trj_a,t,t1); % v(2);

const(5)  = subs(d2trj_a,t,0);
const(6) = subs(d2trj_a,t,t1);

const(7) = subs(d3trj_a,t,0);
const(8) = subs(d3trj_a,t,t1);

const(9) = subs(d4trj_a,t,0);
const(10) = subs(d4trj_a,t,t1);

variables = [A B];

AA = sym(zeros(size(const,2),size(variables,2)));

for i=1:size(const,2)
    for j=1:size(variables,2)
        AA(i,j) = diff(const(i),variables(j));
    end
end

matlabFunction(AA,'file','calc_AA12_','Vars',{[t0 t1]})
