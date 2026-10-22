function task2_filters(I, t, out)
% Task II: Ideal vs Gaussian vs Butterworth, low-pass and high-pass.
[M, N] = size(I);
D = dist_grid(M, N);
F = cfft(I);
types = {'ideal', 'gauss', 'butter'};
names = {'Ideal', 'Gaussian', 'Butterworth n=2'};
D0s = [10 30 60];
kind = {'LPF', 'HPF'};

% 1. Different cutoffs (rows = filter type, columns = cutoff)
fprintf('\n%s  type        D0   PSNR(dB)  edge ratio\n', 'LPF');
for hp = 0:1
    ims = {};  tt = {};
    for a = 1:3
        for b = 1:3
            H = lowpass(types{a}, D, D0s(b), 2);
            if hp, H = 1 - H; end
            G = icfft(F .* H);
            ims{end+1} = G;
            tt{end+1} = sprintf('%s %s D0=%d', names{a}, kind{hp+1}, D0s(b));
            if ~hp
                fprintf('     %-10s %3d   %6.2f    %.3f\n', types{a}, D0s(b), psnr_db(G, I), edge_ratio(G, I));
            end
        end
    end
    show_grid(ims, tt, 3, 3, out, sprintf('%s_T2_filters_%s_cutoffs', t, kind{hp+1}));
end

% 2. Magnitude and phase of original and filtered images (D0 = 30)
for hp = 0:1
    ims = {I, logmag(F), angle(F)};
    tt = {'Original', 'Original magnitude', 'Original phase'};
    for a = 1:3
        H = lowpass(types{a}, D, 30, 2);
        if hp, H = 1 - H; end
        G = cfft(icfft(F .* H));
        ims = [ims, {icfft(G), logmag(G), angle(G)}];
        tt = [tt, {[names{a} ' ' kind{hp+1}], 'Magnitude', 'Phase'}];
    end
    show_grid(ims, tt, 4, 3, out, sprintf('%s_T2_filters_%s_mag_phase', t, kind{hp+1}));
end

% 3. Butterworth order: higher n = sharper cutoff = more ringing
orders = [1 2 4 8 20];
ims = {};  tt = {};
fprintf('\nButterworth LPF D0=30  order  PSNR(dB)  edge ratio  ringing\n');
for hp = 0:1
    for n = orders
        H = lowpass('butter', D, 30, n);
        if hp, H = 1 - H; end
        G = icfft(F .* H);
        ims{end+1} = G;  tt{end+1} = sprintf('Butterworth %s n=%d', kind{hp+1}, n);
        if ~hp
            fprintf('                       %3d    %6.2f    %.3f     %.3f\n', n, ...
                psnr_db(G, I), edge_ratio(G, I), ringing(H));
        end
    end
end
fprintf('Ringing = largest negative side lobe / main peak of the spatial kernel (0 = none)\n');
fprintf('Ideal LPF D0=30 ringing %.3f, Gaussian LPF D0=30 ringing %.3f\n', ...
    ringing(lowpass('ideal', D, 30)), ringing(lowpass('gauss', D, 30)));
show_grid(ims, tt, 2, numel(orders), out, [t '_T2_filters_butter_orders']);

% 4. Transition profile vs spatial kernel (explains ringing)
d = 0:0.5:100;
f = figure('Position', [50 50 1100 420]);
subplot(1, 2, 1); hold on;
plot(d, lowpass('ideal', d, 30), 'k', 'LineWidth', 1.5);
plot(d, lowpass('gauss', d, 30), 'b', 'LineWidth', 1.5);
for n = [1 2 8]
    plot(d, lowpass('butter', d, 30, n), '--', 'LineWidth', 1.2);
end
legend('Ideal', 'Gaussian', 'Butter n=1', 'Butter n=2', 'Butter n=8');
xlabel('D(u,v)'); ylabel('H'); title('Filter transition profile (D0 = 30)');
subplot(1, 2, 2); hold on;
r = floor(M/2) + 1;  c = floor(N/2) + 1;  w = 40;
hs = {lowpass('ideal', D, 30), lowpass('gauss', D, 30), ...
      lowpass('butter', D, 30, 1), lowpass('butter', D, 30, 2), lowpass('butter', D, 30, 8)};
for k = 1:numel(hs)
    h = fftshift(real(ifft2(ifftshift(hs{k}))));
    p = h(r, c-w:c+w);
    plot(-w:w, p / max(p), 'LineWidth', 1.2);
end
legend('Ideal', 'Gaussian', 'Butter n=1', 'Butter n=2', 'Butter n=8');
xlabel('x'); ylabel('h(x) normalized'); title('Spatial kernel (center row): side lobes = ringing');
save_fig(f, out, [t '_T2_filters_profiles']);

% 5. Noise: LPF smooths it, HPF amplifies it relative to the signal
sigma = 0.1;
In = I + sigma * randn(M, N);
Fn = cfft(In);
ims = {In};  tt = {sprintf('Noisy (sigma=%.2f), PSNR %.1f', sigma, psnr_db(In, I))};
for a = 1:3
    G = icfft(Fn .* lowpass(types{a}, D, 50, 2));
    ims{end+1} = G;  tt{end+1} = sprintf('%s LPF D0=50, PSNR %.1f', names{a}, psnr_db(G, I));
end
fprintf('\nHPF D0=30 on noisy image  type   SNR in (dB)  SNR out (dB)\n');
snr_in = 10 * log10(var(I(:)) / sigma^2);
for a = 1:3
    H = 1 - lowpass(types{a}, D, 30, 2);
    Gc = icfft(F .* H);
    Gn = icfft(Fn .* H);
    snr_out = 10 * log10(var(Gc(:)) / var(Gn(:) - Gc(:)));
    fprintf('                          %-6s   %6.2f      %6.2f\n', types{a}, snr_in, snr_out);
    ims{end+1} = Gn;  tt{end+1} = sprintf('%s HPF noisy, SNR %.1f dB', names{a}, snr_out);
end
ims{end+1} = icfft(F .* (1 - lowpass('gauss', D, 30)));  tt{end+1} = 'Gaussian HPF clean';
show_grid(ims, tt, 2, 4, out, [t '_T2_filters_noise']);
end

function r = ringing(H)
% Size of the largest negative side lobe of the spatial kernel, relative to its peak.
h = real(ifft2(ifftshift(H)));
r = max(0, -min(h(:))) / max(h(:));
end
