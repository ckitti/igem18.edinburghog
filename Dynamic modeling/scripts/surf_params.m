clear

%% Initial species concentrations (uM)
S1 = 0; S3 = 0; S5 = 0; S7 = 0;
S2 = 0; S4 = 0; S6 = 0; S8 = 0;
Pr = 1E3;
Ac = 1E3;
CoA = 1000;
NADPH = 1000; NADP = 0;

c = 1;
e0 = c; e1 = c; e1b = c;
e2 = c; e3 = c;
time = 12 *3600; %run time
bins = 30;
mat1 = zeros(bins+1);
mat2 = zeros(bins+1);

%% Compute solutions
for i=1:length(mat1)
    for j=1:length(mat1)
        
        % select species/parameters to vary
        %Ac = 50*(i);
        % Pr = 35*(i-1);
        %CoA = 50*(j);
        % e1 = 0.01*(i-1);
        % e1b = 0.01*(j-1);
        k1 = 10^((i-1)/5 - 3);
        k2 = 10^((j-1)/5 - 3);
        
        %Ac = 10^((i-1)/10);
        %Pr = 10^((j-1)/10);
        
        % select function to run
        %[t1, x1] = ode15s(@fPHBdt, 0:1:time, [Ac CoA S1 S3 S5 S7]);
        %[t2, x2] = ode15s(@fPHBdt_PhaA_rr, 0:1:time, [Ac CoA S1 S3 S5 S7]);
        %[t3, x3] = ode15s(@fPHBdt_rr, 0:1:time, [Ac CoA NADPH NADP S1 S3 S5 S7]);
        %[t4, x4] = ode15s(@fPHBVdt_rr, 0:1:time, [Ac Pr CoA S1 S2 S3 S4 S5 S6 S7 S8]);
        [t4, x4] = ode15s(@(t,x)fPHBVdt_rr_k(t,x,k1,k2), 0:1:time, [Ac Pr CoA S1 S2 S3 S4 S5 S6 S7 S8]);
        
        %[t4, x4] = ode15s(@fPHBdt_v2, 0:1:time, [Ac CoA NADPH NADP S1 S3 S5 S7]);

        %[t4, x4] = ode15s(@(t,x)fPHBVdt_v2(t,x), 0:1:time, [Ac Pr CoA NADPH NADP S1 S2 S3 S4 S5 S6 S7 S8]);

        
        % compiling results into arrays
        mat1(i,j) = x4(end,end); %PHV        
        mat2(i,j) = x4(end,(end-1)); %PHB
        
        smat = mat1 + mat2; %PHBV total
        nmat = mat1 ./ smat; %PHV/PHBV
        
        %[t4, x4] = ode15s(@(t,x)fPHBVdt_rr_variable_e(t,x,e0, e1, e1b, e2, e3),0:1:time, [Ac Pr CoA S1 S2 S3 S4 S5 S6 S7 S8] );

    end
    i
end

figure;
subplot(2,2,1);
surf(mat1);
subplot(2,2,2);
surf(mat2);
subplot(2,2,3);
surf(smat);
subplot(2,2,4);
surf(nmat)