function show_spectra(ims, titles, out, name)
% For each image: image, log magnitude and phase in one column. Saves results/<name>.png.
n = numel(ims);
cells = {};  tt = {};
for k = 1:n
    cells{k} = ims{k};               tt{k} = titles{k};
    F = cfft(ims{k});
    cells{n + k} = logmag(F);        tt{n + k} = 'Magnitude (log)';
    cells{2*n + k} = angle(F);       tt{2*n + k} = 'Phase';
end
show_grid(cells, tt, 3, n, out, name);
end
