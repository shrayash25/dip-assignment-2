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
r = floor(M/2) + 1;
plot(logmag(F(r, :)), 'k'); hold on;
plot(logmag(Fc(r, :)), 'r'); plot(logmag(Fe(r, :)), 'b');
legend('Original', 'Compressed', 'Expanded'); xlabel('v'); ylabel('log(1+|F|)');
title('Center row of spectrum');
save_fig(f, out, [t '_T1_4_scaling_profile']);

% 5. Rotation: smooth circular window removes border effects, so only content rotates
D = dist_grid(M, N);  R = min(M, N) / 2;
W = 0.5 * (1 + cos(pi * min(D / R, 1)));
Iw = I .* W;
angs = [0 30 45 60 90];
ims = {};  tt = {};
for a = angs
    Ir = imrotate(Iw, a, 'bilinear', 'crop');
    ims{end+1} = Ir;  tt{end+1} = sprintf('Rotated %d deg', a);
end
for a = 1:numel(angs)
    ims{end+1} = logmag(cfft(ims{a}));  tt{end+1} = sprintf('Spectrum %d deg', angs(a));
end
show_grid(ims, tt, 2, numel(angs), out, [t '_T1_5_rotation']);
end

function Y = fit_center(X, M, N)
% Pad with zeros or crop so X sits centered in an M x N canvas.
Y = zeros(M, N);
[h, w] = size(X);
r = min(h, M);  c = min(w, N);
xr = floor((h - r)/2);  xc = floor((w - c)/2);
yr = floor((M - r)/2);  yc = floor((N - c)/2);
Y(yr+1:yr+r, yc+1:yc+c) = X(xr+1:xr+r, xc+1:xc+c);
end

function r = energy_radius(F)
% Radius that holds 90% of the non-DC spectral energy.
P = abs(F).^2;
D = dist_grid(size(F, 1), size(F, 2));
P(D == 0) = 0;
[d, idx] = sort(D(:));
e = cumsum(P(idx)) / sum(P(:));
r = d(find(e >= 0.9, 1));
end
