function task2_bands(I, t, out, cutoffs)
% Task II: split image into low, mid, high bands and rebuild each.
% cutoffs: rows of [c1 c2]; low: D<=c1, mid: c1<D<=c2, high: D>c2.
if nargin < 4, cutoffs = [10 40; 30 90]; end
[M, N] = size(I);
D = dist_grid(M, N);
F = cfft(I);

for k = 1:size(cutoffs, 1)
    c1 = cutoffs(k, 1);  c2 = cutoffs(k, 2);
    Hl = double(D <= c1);
    Hm = double(D > c1 & D <= c2);
    Hh = double(D > c2);
    L = icfft(F .* Hl);  Mi = icfft(F .* Hm);  Hi = icfft(F .* Hh);
    S = L + Mi + Hi;
    fprintf('Bands [%d %d]: max |low+mid+high - original| = %.2e\n', c1, c2, max(abs(S(:) - I(:))));
    show_grid({I, L, Mi, Hi, S, logmag(F .* Hl), logmag(F .* Hm), logmag(F .* Hh)}, ...
        {'Original', sprintf('Low  D<=%d', c1), sprintf('Mid  %d<D<=%d', c1, c2), ...
         sprintf('High D>%d', c2), 'Low+Mid+High', 'Low band spectrum', ...
         'Mid band spectrum', 'High band spectrum'}, ...
        2, 4, out, sprintf('%s_T2_bands_%d_%d', t, c1, c2));
end
end
