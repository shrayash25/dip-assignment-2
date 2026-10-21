function show_grid(ims, titles, r, c, out, name)
% Show images in an r x c grid and save the figure as results/<name>.png.
f = figure('Position', [50 50 280*c 300*r]);
for k = 1:numel(ims)
    subplot(r, c, k);
    imshow(ims{k}, []);
    title(titles{k}, 'Interpreter', 'none', 'FontSize', 8);
end
save_fig(f, out, name);
end
