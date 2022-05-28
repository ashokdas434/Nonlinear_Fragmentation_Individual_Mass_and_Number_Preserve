function [x,R,del_x, L] = Lin_Grids(x_min, x_max, I, grid_mul)

R = linspace(x_min,x_max,I+1);

%% Division into smaller cells using multiplier
if grid_mul >1
    R1 = R(1);
    for i=1:I
        del = R(i+1)-R(i);
        for j=1:grid_mul
            R1 = [R1, R1(end)+del/grid_mul];
        end
    end
    R = R1;
    L = length(R)-1;  % No. of cells
else
    L = I;
end
%%
R(L+1) = R(L+1) + (R(L+1)-R(L));

%%
x = zeros(1,L);
for i=1:L
    x(i) = 0.5*(R(i)+R(i+1));
end
del_x = (x_max - x_min)/L * ones(1,L);
end