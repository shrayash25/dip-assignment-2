function D = dist_grid(M, N)
% Distance of each frequency sample from the center of a centered spectrum.
[v, u] = meshgrid((1:N) - floor(N/2) - 1, (1:M) - floor(M/2) - 1);
D = sqrt(u.^2 + v.^2);
end
