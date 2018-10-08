clear

%% Initial species concentration (uM)
S1 = 0; S3 = 0; S5 = 0; S7 = 0;
S2 = 0; S4 = 0; S6 = 0; S8 = 0;
Pr = 1E3;
Ac = 1E3;
CoA = 1000;
NADPH = 1000; NADP = 0;

c = 1;
e0 = c; e1 = c; e1b = c;
e2 = c; e3 = c;
time = 3600*12; %run time
bins = 10;
mat1 = zeros(bins);
mat2 = zeros(bins);

%% Compute solutions
for i=1:length(mat1)
    for j=1:length(mat1)
        %a(i,j) = sin(2*pi*i/100) - cos(2*pi*j/100);
        % Ac = 100 * i;
        Pr = 100 * j ;
        CoA = 100 * i;
        %e1 = 0.01*(i-1);
        %e1b = 0.01*(j-1);
        
        %[t1, x1] = ode15s(@fPHBdt, [0 time], [Ac CoA S1 S3 S5 S7]);
        %[t2, x2] = ode15s(@fPHBdt_PhaA_rr, [0 time], [Ac CoA S1 S3 S5 S7]);
        %[t3, x3] = ode15s(@fPHBdt_rr, [0 time], [Ac CoA NADPH NADP S1 S3 S5 S7]);
        %[t4, x4] = ode15s(@fPHBVdt_rr, 0:1:time, [Ac Pr CoA 0 0 0 0 0 0 0 0]);
        %[t4, x4] = ode15s(@(t,x)fPHBVdt_rr_variable_e(t,x,e0, e1, e1b, e2, e3),0:1:time, [Ac Pr CoA 0 0 0 0 0 0 0 0] );
        %[t4, x4] = ode15s(@fPHBVdt_rr_v2, 0:1:time, [Ac Pr CoA e0 e1 e1b e2 e3 0 0 0 0 0 0 0 0]);
        [t4, x4] = ode15s(@(t,x)fPHBVdt_v2(t,x), 0:1:time, [Ac Pr CoA NADPH NADP S1 S2 S3 S4 S5 S6 S7 S8]);

        mat1(i,j) = x4(end,end);        
        mat2(i,j) = x4(end,(end-1));
        
        smat = mat1 + mat2;
        nmat = mat1 ./ smat;
    end
    i
end
figure;
subplot(4,1,1);
surf(mat1)
subplot(4,1,2);
surf(mat2)
subplot(4,1,3);
surf(smat)
subplot(4,1,4);
surf(nmat)
