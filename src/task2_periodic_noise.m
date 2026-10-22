function task2_periodic_noise(I, t, A, B, out)
% Task II: add periodic interference, find its peaks, remove with notch and band-reject.
[M, N] = size(I);
D = dist_grid(M, N);
[x, y] = meshgrid(0:N-1, 0:M-1);
noise = @(u0, v0, amp) amp * sin(2*pi*(u0*x/N + v0*y/M));

% 1. Different interference frequencies (u0, v0) in cycles per image
freqs = [40 0; 30 30; A*10+5 B*10+25];
amp = 0.2;
ims = {};  tt = {};
fprintf('\nPeriodic noise  (u0,v0)  amp   PSNR noisy  notch  band-reject\n');
for k = 1:size(freqs, 1)
    In = I + noise(freqs(k, 1), freqs(k, 2), amp);
    [Gn, Gb, pk] = denoise(In, D, 5, 8);
    fprintf('                (%d,%d)  %.2f   %6.2f   %6.2f   %6.2f\n', freqs(k, :), amp, ...
        psnr_db(In, I), psnr_db(Gn, I), psnr_db(Gb, I));
    fprintf('                  peaks found at (row,col) offsets (%d,%d) and (%d,%d)\n', pk(1, :), pk(2, :));
    if k == 1
        show_spectra({I, In, Gn, Gb}, {'Original', 'Noisy', 'Notch reject', 'Band reject'}, ...
            out, [t '_T2_noise_spectra']);
    end
    Fn = logmag(cfft(In));
    ims = [ims, {In, mark_peaks(Fn, pk), Gn, Gb}];
    tt = [tt, {sprintf('Noise (%d,%d) amp %.2f', freqs(k, :), amp), 'Spectrum, detected peaks', ...
               sprintf('Notch reject, PSNR %.1f', psnr_db(Gn, I)), ...
               sprintf('Band reject, PSNR %.1f', psnr_db(Gb, I))}];
end
show_grid(ims, tt, 3, 4, out, [t '_T2_noise_frequencies']);

% 2. Different amplitudes
amps = [0.05 0.2 0.5];
ims = {};  tt = {};
for a = amps
    In = I + noise(40, 0, a);
    Gn = denoise(In, D, 5, 8);
    ims = [ims, {In}];  tt = [tt, {sprintf('amp %.2f, PSNR %.1f', a, psnr_db(In, I))}];
    ims = [ims, {Gn}];  tt = [tt, {sprintf('Notch, PSNR %.1f', psnr_db(Gn, I))}];
end
show_grid(ims, tt, 3, 2, out, [t '_T2_noise_amplitudes']);

% 3. Different filter bandwidths: narrow may miss, wide removes image detail
In = I + noise(30, 30, 0.2);
ims = {};  tt = {};
for w = [2 5 15]
    Gn = denoise(In, D, w, 8);
    ims{end+1} = Gn;  tt{end+1} = sprintf('Notch radius %d, PSNR %.1f', w, psnr_db(Gn, I));
end
for w = [3 10 30]
    [~, Gb] = denoise(In, D, 5, w);
    ims{end+1} = Gb;  tt{end+1} = sprintf('Band-reject width %d, PSNR %.1f', w, psnr_db(Gb, I));
end
show_grid(ims, tt, 2, 3, out, [t '_T2_noise_bandwidths']);
end

function [Gn, Gb, pk] = denoise(In, D, r, W)
% Find the two strongest off-center peaks, then apply notch and band-reject filters.
F = cfft(In);
[M, N] = size(In);
pk = find_peaks(F, D, 2);
c = [floor(M/2) + 1, floor(N/2) + 1];
[V, U] = meshgrid(1:N, 1:M);
H = ones(M, N);
for k = 1:size(pk, 1)
    Dk = sqrt((U - c(1) - pk(k, 1)).^2 + (V - c(2) - pk(k, 2)).^2);
    Dm = sqrt((U - c(1) + pk(k, 1)).^2 + (V - c(2) + pk(k, 2)).^2);
    H = H .* (1 ./ (1 + (r ./ (Dk + eps)).^4)) .* (1 ./ (1 + (r ./ (Dm + eps)).^4));
end
Gn = icfft(F .* H);
R0 = norm(pk(1, :));
Gb = icfft(F .* bandreject('butter', D, R0, W, 2));
end

function pk = find_peaks(F, D, k)
% Offsets (row, col) from center of the k sharpest peaks outside the low-frequency area.
% Subtracting a local median keeps isolated spikes and removes the smooth image spectrum.
L = abs(F) - medfilt2(abs(F), [9 9]);
L(D < 15) = 0;
[M, N] = size(F);
c = [floor(M/2) + 1, floor(N/2) + 1];
pk = zeros(k, 2);
for i = 1:k
    [~, idx] = max(L(:));
    [r, cc] = ind2sub([M N], idx);
    pk(i, :) = [r cc] - c;
    L(max(r-3, 1):min(r+3, M), max(cc-3, 1):min(cc+3, N)) = 0;
end
end

function X = mark_peaks(L, pk)
% Draw small boxes around detected peaks for display.
X = L / max(L(:));
[M, N] = size(X);
c = [floor(M/2) + 1, floor(N/2) + 1];
for i = 1:size(pk, 1)
    r = c(1) + pk(i, 1);  cc = c(2) + pk(i, 2);
    rr = max(r-6, 1):min(r+6, M);  cr = max(cc-6, 1):min(cc+6, N);
    X(rr, [cr(1) cr(end)]) = 1;
    X([rr(1) rr(end)], cr) = 1;
end
end
