%%
clear all; clc
%% Parameters
nx = 16; ny = 32;
[rx,ry] = meshgrid(linspace(0, nx - 1, nx), linspace(0, ny - 1, ny));
r = sqrt(rx.^2 + ry.^2);
dx = 0.1; dy = 0.01;

thick = (1 : 1 : 20);

scanning_nx = 5; scanning_ny = 7;
scanning_step = 1;
scanning_x0 = 1; scanning_y0 = 10;
scanning_periodic = false;
xs = (scanning_x0 : scanning_step : scanning_x0 + (scanning_nx - 1) * scanning_step);
ys = (scanning_y0 : scanning_step : scanning_y0 + (scanning_ny - 1) * scanning_step);
[X,Y] = meshgrid(xs, ys);

x = xs(1); y = ys(1);
R = sqrt(X.^2 + Y.^2);
%% EWRS results
clear results
results.images = zeros(nx, ny, size(xs, 2), size(ys, 2), size(thick, 2));
results.thick = thick;
results.thicknesses = {};
results.xs = xs;
results.ys = ys;
results.dx = dx;
results.dy = dy;

for i = 1:size(xs,2)
    x = xs(i);
    for j = 1:size(ys,2)
        y = ys(j);
        results.thicknesses{i, j} = thick;
        for t = 1:size(thick, 2)
            %fprintf('i=%d, j=%d:\n\tSize before assignment: (%d, %d, %d, %d, %d)\n', i, j, size(results.images))
            results.images(:, :, i, j, t) = transpose(r * thick(t) + x*y);
            %fprintf('\tSize after assignment: (%d, %d, %d, %d, %d)\n', size(results.images))
        end
    end
end

save("ewrs_test_results.ecmat", 'results', '-v7.3');

%% SCBED
clear results
results.x=x;
results.y=y;
results.images = zeros(nx, ny, size(thick, 2));
results.thick = thick;
results.dx = dx;
results.dy = dy;
   
for t = 1:size(thick, 2)
    results.images(:, :, t) = transpose(r * thick(t));
end
save("cbed_test_results.ecmat", 'results', '-v7.3');
%% SCBED
clear results
results.xs=xs;
results.ys=ys;
results.images = zeros(nx, ny, size(xs, 2), size(ys, 2), size(thick, 2));
results.thick = thick;
results.thicknesses = {};
results.dx = dx;
results.dy = dy;

for i = 1:size(xs,2)
    x = xs(i);
    for j = 1:size(ys,2)
        y = ys(j);
        results.thicknesses{i, j} = thick;
        for t = 1:size(thick, 2)
            results.images(:, :, i, j, t) = transpose(r * thick(t) + x*y);
        end
    end
end
save("scbed_test_results.ecmat", 'results', '-v7.3');
%% HRTEM
clear results

results.images = zeros(nx, ny, size(thick,2));
results.thick = thick;

results.dx = dx;
results.dy = dy;
   
for t = 1:size(thick, 2)
    results.images(:, :, t) = transpose(r * thick(t));
end
save("hrtem_test_results.ecmat", 'results', '-v7.3');
%% STEM
clear results

results.scanning_periodic = scanning_periodic;
results.scanning_ns = scanning_nx;
results.scanning_x0 = xs(1);
results.scanning_y0 = ys(1);
if scanning_periodic
    results.scanning_xe = xs(length(xs)-1);  
    results.scanning_ye = ys(length(xs)-1);
    results.images = zeros(size(thick,2), scanning_nx-1, scanning_ny-1);
else
    results.scanning_xe = xs(length(xs));
    results.scanning_xe = xs(length(xs));
    results.images = zeros(size(thick,2), scanning_nx, scanning_ny);
end    
results.thick = thick;
results.dx = dx;
results.dy = dy;
   
for t = 1:size(thick, 2)
    results.images(t, :, :) = transpose(R * thick(t));
end

save("stem_test_results.ecmat", 'results', '-v7.3');
