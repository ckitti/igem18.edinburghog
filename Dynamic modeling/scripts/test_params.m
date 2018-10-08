clear
S1 = 0; S3 = 0; S5 = 0; S7 = 0;
S2 = 0; S4 = 0; S6 = 0; S8 = 0;
Pr = 1000;
Ac = 1000;
CoA = 1000;
e0 = 0.1; e1 = 0.1; e1b = 0.1; e2 = 0.1; e3 = 0.1;
time = 720;
mat1 = zeros(11);
%mat2 = zeros(101);
 

for i=1:length(mat1)
    for j=1:length(mat1)
        for h=1:length(mat1)
            Ac = 250*(i-1);
            Pr = 250*(j-1);
            CoA = 250*(h-1);
            %[t1, x1] = ode15s(@fPHBdt, [0 time], [Ac CoA S1 S3 S5 S7]);
            %[t2, x2] = ode15s(@fPHBdt_PhaA_rr, [0 time], [Ac CoA S1 S3 S5 S7]);
            %[t3, x3] = ode15s(@fPHBdt_rr, [0 time], [Ac CoA NADPH NADP S1 S3 S5 S7]);
            [t4, x4] = ode15s(@fPHBVdt_rr, 0:1:time, [Ac Pr CoA 0 0 0 0 0 0 0 0]);
            %[t4, x4] = ode15s(@(t,x)fPHBVdt_rr_variable_e(t,x,e0, e1, e1b, e2, e3),0:1:time, [Ac Pr CoA 0 0 0 0 0 0 0 0] );
            dat(i,j,h) = x4(end,end);        
            % mat2(i,j) = x4(end,(end-1));
        
            %nmat = mat1 ./ mat2;
        end
        i
    end
end
%figure;
%subplot(2,1,1);
%surf(mat1)
%subplot(2,1,2);
%plot(t2/3600,x2(:,end))
%surf(mat2)