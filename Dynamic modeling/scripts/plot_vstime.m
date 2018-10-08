clear
S1 = 0; S3 = 0; S5 = 0; S7 = 0; 
S2 = 0; S4 = 0; S6 = 0; S8 = 0;
NADPH = 1000; NADP = 0;
time = 72 * 3600;
CoA = 1000;
Ac = 1E3; Pr = 1E3;
e0 = 0.1; e1 = 0.1; e1b = 0.1; e2 = 0.1; e3 = 0.1;

%kcat
%k1 = 26.46;

%specific activity
k1 = 100;

for i=1:5
    Ac = 1000;
    Pr = 500 * i;
    %CoA = 1000;
    %e1b = e1 * 10^(2*i-1)
    %k2 = k1 * 10^(i-3)

    
    %[t1, x1] = ode15s(@fPHBdt, 0:1:time, [Ac CoA S1 S3 S5 S7]);
    %[t2, x2] = ode15s(@fPHBdt_PhaA_rr, 0:1:time,[Ac CoA S1 S3 S5 S7]);
    %[t3, x3] = ode15s(@fPHBdt_rr,0:1:time, [Ac CoA NADPH NADP S1 S3 S5 S7]);
    
    %[t4, x4] = ode15s(@fPHBVdt_rr, 0:1:time, [Ac Pr CoA 0 0 0 0 0 0 0 0]);
    %[t4, x4] = ode15s(@(t,x)fPHBVdt_rr_variable_e(t,x,e0, e1, e1b, e2, e3), 0:1:time,[Ac Pr CoA 0 0 0 0 0 0 0 0] );
%    [t4, x4] = ode15s(@(t,x)fPHBVdt_rr_k(t,x,k1,k2), 0:1:time, [Ac Pr CoA S1 S2 S3 S4 S5 S6 S7 S8]);
    [t4, x4] = ode15s(@(t,x)fPHBVdt_rr(t,x), 0:1:time, [Ac Pr CoA S1 S2 S3 S4 S5 S6 S7 S8]);


    %[t4, x4] = ode15s(@(t,x)fPHBdt_v2(t,x), 0:1:time, [Ac CoA NADPH NADP S1 S3 S5 S7]);
    
%    [t4, x4] = ode15s(@(t,x)fPHBVdt_v2(t,x), 0:1:time, [Ac Pr CoA NADPH NADP S1 S2 S3 S4 S5 S6 S7 S8]);
%    [t4, x4] = ode15s(@(t,x)fPHBVdt_v2(t,x,k1,k2), 0:1:time, [Ac Pr CoA NADPH NADP S1 S2 S3 S4 S5 S6 S7 S8]);

    %y1(:,i) = x1(:,end);
    %y2(:,i) = x2(:,end);
    %y3(:,i) = x3(:,end);
    %y4(:,i) = x4(:,end);
    i
    y1(:,i) = x4(:,end);
    y2(:,i) = x4(:,(end-1));
    %y(:,i) = x4(:,end)./x4(:,(end-1));
    sy(:,i) = x4(:,end) + x4(:, (end-1));
    ny(:,i) = x4(:,end)./(x4(:,end)+x4(:,(end-1)));
    [X, I] = max(ny(:,i)) % max height and the corresponding index (time in s)
    
    
    %s8(:,i) = x4(:,end);
    %s7(:,i) = x4(:,(end-1));
    %s6(:,i) = x4(:,(end-2));
    %s5(:,i) = x4(:,(end-3));
    %s4(:,i) = x4(:,(end-4));
    %s3(:,i) = x4(:,(end-5));
    %s2(:,i) = x4(:,(end-6));
    %s1(:,i) = x4(:,(end-7));
    %coa(:,i) = x4(:,(end-8));
    %pr(:,i) = x4(:,(end-9));
    %ac(:,i) = x4(:,(end-10));

    %v8 = diff([s8; zeros(1,5)])
    %final
    % v8 = [diff(s8); zeros(1,5)]
end


% fit