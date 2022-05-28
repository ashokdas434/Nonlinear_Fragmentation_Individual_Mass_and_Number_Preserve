x_min = 0; x_max = 1;  % property limits
I     = 250; % Number of intervals
len_T = 11;
T = linspace(0,5,len_T);  % [sec] Process time

[x,R,del_x] = Grids2(x_min, x_max, I);

%% Analytical Solution
f = @(x,t) exp(-x*t).*(2.*t+ t.^2 .*(1-x)); % Analytical solution without the \delta(x-1)*exp(-t) term
xf = @(x,t) exp(-x*t).*(2.*t.*x+ t.^2 .*(x-x.^2));

N_ana = zeros(len_T,I); M_dist = zeros(len_T,I); % Initialize

for j=1:len_T
    for i=1:I-1
        N_ana(j,i) = integral(@(x) f(x,T(j)),R(i),R(i+1));
        M_dist(j,i) = integral(@(x) xf(x,T(j)),R(i),R(i+1));
    end
    N_ana(j,I) = integral(@(x) f(x,T(j)),R(I),x_max) + exp(-T(j));
    M_dist(j,I) = integral(@(x) xf(x,T(j)),R(I),x_max) + exp(-T(j));
end
 
N_tot_ana = sum(N_ana,2); % Total particles
 
M_tot = sum(M_dist,2); % Total mass

%% Test
M_tot_test = N_ana*x';