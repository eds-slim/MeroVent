fid = fopen('./../R/mic2.dat');
d = textscan(fid, '%d%f%d', 'Delimiter', ',', 'HeaderLines', 1);
fclose(fid)

dd = cell2table(d);

heatmap(dd, 'd1', 'd2', 'ColorVariable', 'd3')

%%

surf(d{1},d{2},d{3})

%%
x=d{1};
y=d{2};
z=d{3};