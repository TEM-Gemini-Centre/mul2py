function [results_struct] = create_dummy_results()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
results_struct = struct();
nx = 128;
ny = 128;
t = 10;
nX = 3;
nY = 4;

results_struct.images = rand(nx, ny, t, nX, nY);
%results_struct.images = rand(nx, ny, t);
results_struct.title = 'test';
results_struct.xs = [1:nx];
results_struct.ys = [1:nY];
input_multem = HRTEM_setup('C:\Users\emilc\OneDrive - NTNU\NORTEM\Projects\mul2py\mul2py\examples\Models\Al_10x10x20.mat', 'instrument', '2100F');
results_struct.input = input_multem.toStruct();
results_struct.dz = 2.025;
results_struct.dx = 1.1;
results_struct.dy = 1.2;
results_struct.thick = [results_struct.dz:results_struct.dz:results_struct.dz*size(results_struct.images, 3)];
results_struct.axes = setup_axes(results_struct);
end

