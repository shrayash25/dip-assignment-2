function save_fig(f, out, name)
% Save figure as PNG and close it.
saveas(f, fullfile(out, [name '.png']));
close(f);
end
