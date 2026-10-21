function r = edge_ratio(X, ref)
% Mean gradient magnitude of X divided by that of ref (1 = edges kept).
[gx, gy] = gradient(X);
[rx, ry] = gradient(ref);
r = mean(sqrt(gx(:).^2 + gy(:).^2)) / mean(sqrt(rx(:).^2 + ry(:).^2));
end
