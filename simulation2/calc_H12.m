function H_ = calc_H12(in1)
%CALC_H12
%    H_ = CALC_H12(IN1)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    31-Aug-2017 16:54:44

t1 = in1(:,2);
t3 = t1.^2;
t4 = t3.^2;
t5 = t4.^2;
t6 = t3.*1.44e3;
t7 = t1.*t3.*2.88e3;
t8 = t4.*1.08e4;
t9 = t4.*5.04e3;
t10 = t1.*t4.*2.016e4;
t11 = t3.*t4.*5.04e4;
t12 = t1.*t4.*8.064e3;
t13 = t3.*t4.*3.36e4;
t14 = t1.*t3.*t4.*8.64e4;
t15 = t5.*1.764e5;
t16 = t3.*t4.*1.2096e4;
t17 = t1.*t3.*t4.*5.184e4;
t18 = t5.*1.3608e5;
t19 = t1.*t5.*2.8224e5;
t20 = t3.*t5.*5.08032e5;
t21 = t1.*t3.*t4.*1.728e4;
t22 = t5.*7.56e4;
t23 = t1.*t5.*2.016e5;
t24 = t3.*t5.*4.2336e5;
t25 = t1.*t3.*t5.*7.697454545454545e5;
t26 = t4.*t5.*1.27008e6;
t27 = t5.*2.376e4;
t28 = t1.*t5.*1.056e5;
t29 = t3.*t5.*2.8512e5;
t30 = t1.*t3.*t5.*6.048e5;
t31 = t4.*t5.*1.1088e6;
t32 = t1.*t4.*t5.*1.842313846153846e6;
t33 = t3.*t4.*t5.*2.8512e6;
t34 = t1.*t5.*3.168e4;
t35 = t3.*t5.*1.4256e5;
t36 = t1.*t3.*t5.*3.888e5;
t37 = t4.*t5.*8.316e5;
t38 = t1.*t4.*t5.*1.535261538461538e6;
t39 = t3.*t4.*t5.*2.56608e6;
t40 = t1.*t3.*t4.*t5.*3.99168e6;
t41 = t5.^2;
t42 = t41.*5.8806e6;
H_ = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t1.*5.76e2,t6,t7,t9,t12,t16,t21,t27,t34,0.0,0.0,0.0,0.0,t6,t1.*t3.*4.8e3,t8,t10,t13,t17,t22,t28,t35,0.0,0.0,0.0,0.0,t7,t8,t1.*t4.*2.592e4,t11,t14,t18,t23,t29,t36,0.0,0.0,0.0,0.0,t9,t10,t11,t1.*t3.*t4.*1.008e5,t15,t19,t24,t30,t37,0.0,0.0,0.0,0.0,t12,t13,t14,t15,t1.*t5.*3.136e5,t20,t25,t31,t38,0.0,0.0,0.0,0.0,t16,t17,t18,t19,t20,t1.*t3.*t5.*8.313250909090909e5,t26,t32,t39,0.0,0.0,0.0,0.0,t21,t22,t23,t24,t25,t26,t1.*t4.*t5.*1.953969230769231e6,t33,t40,0.0,0.0,0.0,0.0,t27,t28,t29,t30,t31,t32,t33,t1.*t3.*t4.*t5.*4.18176e6,t42,0.0,0.0,0.0,0.0,t34,t35,t36,t37,t38,t39,t40,t42,t1.*t41.*8.302023529411765e6],[13,13]);
