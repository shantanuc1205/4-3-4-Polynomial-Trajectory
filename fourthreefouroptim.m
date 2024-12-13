function [QP, p_Vals, t_Vals, dp_Vals, ddp_Vals] = fourthreefouroptim(waypoints, times)
    A_eq = zeros(14,14);
    Q = zeros(14);
    b_eq = zeros(14,1);
    p_Vals = [];
    t_Vals = [];
    dp_Vals = [];
    ddp_Vals = [];

    for i = 1:3:(length(waypoints)-3)
        
        P0 = waypoints(i);
        P1 = waypoints(i+1);
        P2 = waypoints(i+2);
        P3 = waypoints(i+3);

        T0 = sum(times(1:i));
        T1 = sum(times(1:i+1));
        T2 = sum(times(1:i+2));
        T3 = sum(times(1:i+3));
        
        t1 = T1 - T0;
        t2 = T2 - T1;
        t3 = T3 - T2;

        %Hessian Setup
        Q(4,4) = 72;
        Q(5,5) = 1152*(t1^2);
        Q(4,5) = 288*t1;
        Q(5,4) = Q(4,5);
        Q(9,9) = 72;
        Q(12,12) = 72;
        Q(13,13) = 1152*(t3^2);
        Q(12,13) = 288*t3;
        Q(13,12) = Q(12,13);

        %A_eq
        A_eq(1,1) = 1; %pos at t = 0
        A_eq(2,2) = 1; %vel at t =0 
        A_eq(3,3) = 2; %acc at t = 0
        A_eq(4,1:5) = [1, t1, t1^2, t1^3, t1^4]; %pos at t1
        % A_eq(4,2:5) = [1, 2*t1, 3*t1^2, 4*t1^3] %final velocity
        % A_eq(5,3:5) = [2, 6*t1, 12*t1^2] %final accel
        A_eq(5,1:6) = [1, t1, t1^2, t1^3, t1^4, -1]; %pos cont
        A_eq(6,2:7) = [1, 2*t1, 3*t1^2, 4*t1^3, 0, -1 ]; %vel cont
        A_eq(7,3:8) = [2, 6*t1, 12*t1^2, 0, 0, -2]; %accel cont
        A_eq(8,6:9) = [1, t2, t2^2, t2^3]; %end pos t2
        A_eq(9,6:10) = [1, t2, t2^2, t2^3, -1]; %pos cont
        A_eq(10,7:11) = [1, 2*t2, 3*t2^2, 0, -1]; %end vel cont 
        A_eq(11,8:12) = [2, 6*t2, 0, 0, -2]; %accel cont
        A_eq(12,10:14) = [1, t3, t3^2, t3^3, t3^4]; %pos end t3
        A_eq(13,11:14) = [1, 2*t3, 3*t3^2, 4*t3^3]; %vel end t3
        A_eq(14,12:14) = [2, 6*t3, 12*t3^2]; %accel end t3

        %b_eq
        b_eq(1) = waypoints(i);
        b_eq(2) = 0;
        b_eq(3) = 0;
        b_eq(4) = waypoints(i+1);
        b_eq(5) = 0;
        b_eq(6) = 0;
        b_eq(7) = 0;
        b_eq(8) = waypoints(i+2);
        b_eq(9) = 0;
        b_eq(10) = 0;
        b_eq(11) = 0;
        b_eq(12) = waypoints(i+3);
        b_eq(13) = 0;
        b_eq(14) = 0;
        A_eq; b_eq;
        QP = quadprog(Q,[],[],[],A_eq,b_eq);
        A_eq = zeros(14,14);
        Q = zeros(14);
        b_eq = zeros(14,1);

        %polynomials
        [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13] = deal(QP(1), QP(2), QP(3), QP(4), QP(5), QP(6), QP(7), QP(8), QP(9), QP(10), QP(11), QP(12), QP(13), QP(14));
        a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13;
        
        t = linspace(T0,T1,30);
        poly1 = a0 + a1*(t-T0) + a2*(t-T0).^2 + a3*(t-T0).^3 + a4*(t-T0).^4;
        dpoly1 = a1 + 2*a2*(t-T0) + 3*a3*(t-T0).^2 + 4*a4*(t-T0).^3;
        ddpoly1 = 2*a2 + 6*a3*(t-T0) + 12*a4*(t-T0).^2;
        p_Vals = [p_Vals,poly1];
        dp_Vals = [dp_Vals, dpoly1];
        ddp_Vals = [ddp_Vals, ddpoly1];
        t_Vals = [t_Vals, t];
%         plot(t,poly1);
%         hold on;

        t = linspace(T1,T2,30);
        poly2 = a5 + a6*(t-T1) + a7*(t-T1).^2 + a8*(t-T1).^3;
        dpoly2 = a6 + 2*a7*(t-T1) + 3*a8*(t-T1).^2;
        ddpoly2 = 2*a7 + 6*a8*(t-T1);
        p_Vals = [p_Vals,poly2];
        dp_Vals = [dp_Vals, dpoly2];
        ddp_Vals = [ddp_Vals, ddpoly2];
        t_Vals = [t_Vals, t];
%         plot(t,poly2);

        t = linspace(T2,T3,30);
        poly3 = a9 + a10*(t-T2) + a11*(t-T2).^2 + a12*(t-T2).^3 + a13*(t-T2).^4;
        dpoly3 = a10 + 2*a11*(t-T2) + 3*a12*(t-T2).^2 + 4*a13*(t-T2).^3;
        ddpoly3 = 2*a11 + 6*a12*(t-T0) + 12*a13*(t-T0).^2;
        p_Vals = [p_Vals,poly3];
        dp_Vals = [dp_Vals, dpoly3];
        ddp_Vals = [ddp_Vals, ddpoly3];
        t_Vals = [t_Vals, t];
%         plot(t,poly3);
%         hold on


    end

end