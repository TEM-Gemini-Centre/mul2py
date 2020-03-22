%% EWRS results
clear all
thick=(2:1:20);
nx=16; ny=16;
xs = (3:1:6);
ys = (4:1:7);

dx = 1/nx;
dy = 1/ny;

results.images = zeros(ny, nx, size(ys, 2), size(xs, 2), size(thick, 2));
results.thick = thick;
results.thicknesses = {};

results.xs = xs;
results.ys = ys;

for i = 1:size(results.images, 4)
    for j = 1:size(results.images, 5)
        results.thicknesses{i, j} = thick;
        for t = 1:size(thick, 2)
            results.images(:, :, j, i, t) = rand(ny, nx);
        end
    end
end
results.dx = dx;
results.dy = dy;

save("ewrs_test_results.ecmat", 'results', '-v7.3');

%% CBED
clear all
thick=(2:1:20);
nx=16; ny=16;
x = 3;
y = 4;

dx = 1/nx;
dy = 1/ny;

results.x=x;
results.y=y;

results.images = zeros(ny, nx, size(thick, 2));
   
for t = 1:size(thick, 2)
    results.images(:, :, t) = rand(ny, nx);
end

results.thick = thick;

results.dx = dx;
results.dy = dy;

save("cbed_test_results.ecmat", 'results', '-v7.3');
%% SCBED
clear all
thick=(2:1:20);
nx=16; ny=16;
xs = (3:1:6);
ys = (4:1:7);

dx = 1/nx;
dy = 1/ny;

results.xs=xs;
results.ys=ys;


results.images = zeros(ny, nx, size(ys, 2), size(xs, 2), size(thick, 2));
results.thick = thick;
results.thicknesses = {};

for i = 1:size(results.images, 3)
    for j = 1:size(results.images, 4)
        x = xs(j);
        y = ys(i);
        results.thicknesses{i, j} = thick;
        for t = 1:size(thick, 2)
            results.images(:, :, j, i, t) = rand(ny, nx);
        end
    end
end
results.dx = dx;
results.dy = dy;

save("scbed_test_results.ecmat", 'results', '-v7.3');
%% HRTEM
clear all
thick=(2:1:20);
nx=16; ny=16;

dx = 1/nx;
dy = 1/ny;

results.images = zeros(ny, nx, size(thick, 2));
   
for t = 1:size(thick, 2)
    results.images(:, :, t) = rand(ny, nx);
end

results.thick = thick;

results.dx = dx;
results.dy = dy;

save("hrtem_test_results.ecmat", 'results', '-v7.3');
%% STEM
clear all
thick=(2:1:20);
nx=16; ny=16;

ns = 5;
x0 = 0;
y0 = 3;
xe = 5;
ye = 8;

dx = 1/nx;
dy = 1/ny;

results.scanning_ns = ns;
results.scanning_x0 = x0;
results.scanning_xe = xe;
results.scanning_y0 = y0;
results.scanning_ye = ye;

results.images = zeros(ny, nx, size(thick, 2));
   
for t = 1:size(thick, 2)
    results.images(:, :, t) = rand(ny, nx);
end

results.thicknesses = thick;

results.dx = dx;
results.dy = dy;

save("stem_test_results.ecmat", 'results', '-v7.3');
