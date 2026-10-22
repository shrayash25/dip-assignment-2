function task1_properties(I, t, A, B, out)
% Task I: spectrum display, fftshift, shifting, scaling, rotation.
[M, N] = size(I);
F = cfft(I);

% 1. Original, log magnitude, phase
show_grid({I, logmag(F), angle(F)}, ...
    {'Original', 'Log magnitude', 'Phase'}, 1, 3, out, [t '_T1_1_spectrum']);

% 2. Why fftshift: fft2 puts DC at (1,1), so low frequencies sit in the corners
F0 = fft2(I);
[~, i0] = max(abs(F0(:)));  [r0, c0] = ind2sub(size(F0), i0);
[~, i1] = max(abs(F(:)));   [r1, c1] = ind2sub(size(F), i1);
fprintf('fftshift: DC at (%d,%d) without shift, at (%d,%d) with shift\n', r0, c0, r1, c1);

% A centered low-pass mask only works on a centered spectrum
H = lowpass('gauss', dist_grid(M, N), 30);
wrong = real(ifft2(F0 .* H));
right = icfft(F .* H);
% Same centering without fftshift: multiply image by (-1)^(x+y)
[x, y] = meshgrid(0:N-1, 0:M-1);
Fm = fft2(I .* (-1).^(x + y));
show_grid({logmag(F0), logmag(F), logmag(Fm), H, wrong, right}, ...
    {'No fftshift (DC in corners)', 'With fftshift', 'fft2 of I.*(-1)^(x+y)', ...
     'Centered Gaussian LPF', 'LPF on unshifted spectrum', 'LPF on shifted spectrum'}, ...
    2, 3, out, [t '_T1_2_fftshift']);

% 3. Shifting: magnitude unchanged, phase gets a linear ramp
dx = A*10 + 5;  dy = B*10 + 5;
Is = circshift(I, [dy dx]);
Fs = cfft(Is);
fprintf('Shift (%d,%d): max |magnitude change| = %.2e\n', dx, dy, max(abs(abs(Fs(:)) - abs(F(:)))));
show_grid({I, logmag(F), angle(F), Is, logmag(Fs), angle(Fs), ...
           abs(abs(Fs) - abs(F)), angle(Fs .* conj(F))}, ...
    {'Original', 'Original magnitude', 'Original phase', ...
     sprintf('Shifted (%d,%d)', dx, dy), 'Shifted magnitude', 'Shifted phase', ...
     'Magnitude difference', 'Phase difference (linear ramp)'}, ...
    3, 3, out, [t '_T1_3_shift']);

% 4. Scaling: compress by A, expand by B, same canvas size
Ic = fit_center(imresize(I, 1/A), M, N);
Ie = fit_center(imresize(I, B), M, N);
Fc = cfft(Ic);  Fe = cfft(Ie);
fprintf('Scaling: 90%% energy radius  original %d, compressed(1/%d) %d, expanded(x%d) %d\n', ...
    round(energy_radius(F)), A, round(energy_radius(Fc)), B, round(energy_radius(Fe)));
show_grid({I, Ic, Ie, logmag(F), logmag(Fc), logmag(Fe)}, ...
    {'Original', sprintf('Compressed by %d', A), sprintf('Expanded by %d', B), ...
     'Original spectrum', 'Compressed: spectrum spreads', 'Expanded: spectrum shrinks'}, ...
    2, 3, out, [t '_T1_4_scaling']);

f = figure('Position', [50 50 800 400]);
[r0, e0] = energy_curve(F);  [rc, ec] = energy_curve(Fc);  [re, ee] = energy_curve(Fe);
plot(r0, e0, 'k', rc, ec, 'r', re, ee, 'b', 'LineWidth', 1.5); hold on;
plot([0 150], [0.9 0.9], 'k--');
xlim([0 150]); grid on;
legend('Original', sprintf('Compressed by %d', A), sprintf('Expanded by %d', B), '90%', ...
    'Location', 'southeast');
xlabel('Radius D(u,v)'); ylabel('Fraction of energy inside radius');
title('Cumulative spectral energy (DC removed)');
save_fig(f, out, [t '_T1_4_scaling_profile']);

% 5. Rotation: window is flat in the middle and tapers to zero at the circle edge,
%    so the square image border does not add fixed lines to the spectrum
D = dist_grid(M, N);  R = min(M, N) / 2;
taper = min(max((D - 0.7*R) / (0.3*R), 0), 1);
Iw = I .* (0.5 * (1 + cos(pi * taper)));
c = [floor(M/2) + 1, floor(N/2) + 1];  z = 64;
center_crop = @(L) L(c(1)-z:c(1)+z, c(2)-z:c(2)+z);
msk = dist_grid(2*z + 1, 2*z + 1) < z - 2;
L0 = center_crop(logmag(cfft(Iw)));
angs = [0 30 45 60 90];
imgs = {};  specs = {};  tt1 = {};  tt2 = {};
fprintf('Rotation: corr(spectrum, original spectrum rotated) vs corr(spectrum, original spectrum)\n');
for a = angs
    Ir = imrotate(Iw, a, 'bilinear', 'crop');
    L = center_crop(logmag(cfft(Ir)));
    Lref = imrotate(L0, a, 'bilinear', 'crop');
    fprintf('  %2d deg: %.3f vs %.3f\n', a, corr2(L(msk), Lref(msk)), corr2(L(msk), L0(msk)));
    imgs{end+1} = Ir;  specs{end+1} = L;
    tt1{end+1} = sprintf('Rotated %d deg', a);
    tt2{end+1} = sprintf('Spectrum %d deg (center zoom)', a);
end
show_grid([imgs specs], [tt1 tt2], 2, numel(angs), out, [t '_T1_5_rotation']);
end

function Y = fit_center(X, M, N)
% Pad with the mean (no artificial border edges) or crop so X sits centered in M x N.
Y = mean(X(:)) * ones(M, N);
[h, w] = size(X);
r = min(h, M);  c = min(w, N);
xr = floor((h - r)/2);  xc = floor((w - c)/2);
yr = floor((M - r)/2);  yc = floor((N - c)/2);
Y(yr+1:yr+r, yc+1:yc+c) = X(xr+1:xr+r, xc+1:xc+c);
end

function [d, e] = energy_curve(F)
% Fraction of non-DC spectral energy inside each radius.
P = abs(F).^2;
D = dist_grid(size(F, 1), size(F, 2));
P(D == 0) = 0;
[d, idx] = sort(D(:));
e = cumsum(P(idx)) / sum(P(:));
end

function r = energy_radius(F)
% Radius that holds 90% of the non-DC spectral energy.
[d, e] = energy_curve(F);
r = d(find(e >= 0.9, 1));
end
