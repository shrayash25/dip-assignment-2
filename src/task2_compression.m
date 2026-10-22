function task2_compression(I, t, out)
% Task II: keep parts of the spectrum, compare error vs fraction of kept coefficients.
[M, N] = size(I);
D = dist_grid(M, N);
F = cfft(I);
total = M * N;

% 1. Keep only central frequencies (radius R)
Rs = [5 10 20 40 80 160];
frac_c = zeros(size(Rs));  err_c = frac_c;  ims = {};  tt = {};
for k = 1:numel(Rs)
    K = D <= Rs(k);
    G = icfft(F .* K);
    frac_c(k) = nnz(K) / total;
    err_c(k) = psnr_db(G, I);
    ims{end+1} = G;
    tt{end+1} = sprintf('R=%d, kept %.2f%%, PSNR %.1f', Rs(k), 100*frac_c(k), err_c(k));
end
show_grid(ims, tt, 2, 3, out, [t '_T2_compress_central']);
fprintf('\nCompression  central disk R   kept %%   PSNR(dB)\n');
fprintf('             %4d           %7.3f   %6.2f\n', [Rs; 100*frac_c; err_c]);

% 2. Remove selected bands [r1 r2]
bands = [0 0; 5 20; 20 60; 60 120; 120 400];
ims = {};  tt = {};
fprintf('\nCompression  removed band   kept %%   PSNR(dB)\n');
for k = 1:size(bands, 1)
    K = ~(D > bands(k, 1) & D <= bands(k, 2));
    G = icfft(F .* K);
    fr = nnz(K) / total;
    fprintf('             (%3d,%3d]     %6.2f   %6.2f\n', bands(k, :), 100*fr, psnr_db(G, I));
    ims{end+1} = G;
    tt{end+1} = sprintf('Removed (%d,%d], kept %.1f%%, PSNR %.1f', bands(k, :), 100*fr, psnr_db(G, I));
end
show_grid(ims, tt, 2, 3, out, [t '_T2_compress_bands']);

% 3. Keep the largest coefficients (best k-term), compare with central low-pass
fracs = [0.001 0.003 0.01 0.03 0.1 0.3];
err_t = zeros(size(fracs));
mag = sort(abs(F(:)), 'descend');
for k = 1:numel(fracs)
    thr = mag(max(1, round(fracs(k) * total)));
    err_t(k) = psnr_db(icfft(F .* (abs(F) >= thr)), I);
end
fprintf('\nCompression  largest coeffs   kept %%   PSNR(dB)\n');
fprintf('                             %7.3f   %6.2f\n', [100*fracs; err_t]);

f = figure('Position', [50 50 700 450]);
semilogx(100*frac_c, err_c, 'o-', 'LineWidth', 1.5); hold on;
semilogx(100*fracs, err_t, 's-', 'LineWidth', 1.5);
grid on; xlabel('Coefficients kept (%)'); ylabel('PSNR (dB)');
legend('Central disk (low-pass)', 'Largest magnitudes', 'Location', 'southeast');
title('Reconstruction quality vs retained coefficients');
save_fig(f, out, [t '_T2_compress_curve']);
end
