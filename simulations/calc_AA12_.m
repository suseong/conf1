function AA = calc_AA12_(in1)
%CALC_AA12_
%    AA = CALC_AA12_(IN1)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    31-Aug-2017 16:48:54

t1 = in1(:,2);
t2 = in1(:,3);
t4 = t1.^2;
t5 = t4.^2;
t6 = t5.^2;
t7 = t2.^2;
t8 = t7.^2;
t9 = t8.^2;
AA = reshape([1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t1,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t4,0.0,0.0,0.0,t1.*2.0,0.0,0.0,2.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t1.*t4,0.0,0.0,0.0,t4.*3.0,0.0,0.0,0.0,t1.*6.0,0.0,6.0,6.0,0.0,0.0,0.0,0.0,0.0,t5,0.0,0.0,0.0,t1.*t4.*4.0,0.0,0.0,0.0,t4.*1.2e1,0.0,0.0,t1.*2.4e1,0.0,2.4e1,2.4e1,0.0,0.0,t1.*t5,0.0,0.0,0.0,t5.*5.0,0.0,0.0,0.0,t1.*t4.*2.0e1,0.0,0.0,t4.*6.0e1,0.0,0.0,t1.*1.2e2,0.0,0.0,t4.*t5,0.0,0.0,0.0,t1.*t5.*6.0,0.0,0.0,0.0,t5.*3.0e1,0.0,0.0,t1.*t4.*1.2e2,0.0,0.0,t4.*3.6e2,0.0,0.0,t1.*t4.*t5,0.0,0.0,0.0,t4.*t5.*7.0,0.0,0.0,0.0,t1.*t5.*4.2e1,0.0,0.0,t5.*2.1e2,0.0,0.0,t1.*t4.*8.4e2,0.0,0.0,t6,0.0,0.0,0.0,t1.*t4.*t5.*8.0,0.0,0.0,0.0,t4.*t5.*5.6e1,0.0,0.0,t1.*t5.*3.36e2,0.0,0.0,t5.*1.68e3,0.0,0.0,t1.*t6,0.0,0.0,0.0,t6.*9.0,0.0,0.0,0.0,t1.*t4.*t5.*7.2e1,0.0,0.0,t4.*t5.*5.04e2,0.0,0.0,t1.*t5.*3.024e3,0.0,0.0,t4.*t6,0.0,0.0,0.0,t1.*t6.*1.0e1,0.0,0.0,0.0,t6.*9.0e1,0.0,0.0,t1.*t4.*t5.*7.2e2,0.0,0.0,t4.*t5.*5.04e3,0.0,0.0,t1.*t4.*t6,0.0,0.0,0.0,t4.*t6.*1.1e1,0.0,0.0,0.0,t1.*t6.*1.1e2,0.0,0.0,t6.*9.9e2,0.0,0.0,t1.*t4.*t5.*7.92e3,0.0,0.0,t5.*t6,0.0,0.0,0.0,t1.*t4.*t6.*1.2e1,0.0,0.0,0.0,t4.*t6.*1.32e2,0.0,0.0,t1.*t6.*1.32e3,0.0,0.0,t6.*1.188e4,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t7,0.0,0.0,0.0,t2.*2.0,0.0,-2.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2.*t7,0.0,0.0,0.0,t7.*3.0,0.0,0.0,t2.*6.0,0.0,-6.0,6.0,0.0,0.0,0.0,0.0,0.0,0.0,t8,0.0,0.0,0.0,t2.*t7.*4.0,0.0,0.0,t7.*1.2e1,0.0,0.0,t2.*2.4e1,0.0,-2.4e1,2.4e1,0.0,0.0,0.0,t2.*t8,0.0,0.0,0.0,t8.*5.0,0.0,0.0,t2.*t7.*2.0e1,0.0,0.0,t7.*6.0e1,0.0,0.0,t2.*1.2e2,0.0,0.0,0.0,t7.*t8,0.0,0.0,0.0,t2.*t8.*6.0,0.0,0.0,t8.*3.0e1,0.0,0.0,t2.*t7.*1.2e2,0.0,0.0,t7.*3.6e2,0.0,0.0,0.0,t2.*t7.*t8,0.0,0.0,0.0,t7.*t8.*7.0,0.0,0.0,t2.*t8.*4.2e1,0.0,0.0,t8.*2.1e2,0.0,0.0,t2.*t7.*8.4e2,0.0,0.0,0.0,t9,0.0,0.0,0.0,t2.*t7.*t8.*8.0,0.0,0.0,t7.*t8.*5.6e1,0.0,0.0,t2.*t8.*3.36e2,0.0,0.0,t8.*1.68e3,0.0,0.0,0.0,t2.*t9,0.0,0.0,0.0,t9.*9.0,0.0,0.0,t2.*t7.*t8.*7.2e1,0.0,0.0,t7.*t8.*5.04e2,0.0,0.0,t2.*t8.*3.024e3,0.0,0.0,0.0,t7.*t9,0.0,0.0,0.0,t2.*t9.*1.0e1,0.0,0.0,t9.*9.0e1,0.0,0.0,t2.*t7.*t8.*7.2e2,0.0,0.0,t7.*t8.*5.04e3,0.0,0.0,0.0,t2.*t7.*t9,0.0,0.0,0.0,t7.*t9.*1.1e1,0.0,0.0,t2.*t9.*1.1e2,0.0,0.0,t9.*9.9e2,0.0,0.0,t2.*t7.*t8.*7.92e3,0.0,0.0,0.0,t8.*t9,0.0,0.0,0.0,t2.*t7.*t9.*1.2e1,0.0,0.0,t7.*t9.*1.32e2,0.0,0.0,t2.*t9.*1.32e3,0.0,0.0,t9.*1.188e4],[17,26]);
