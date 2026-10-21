function task2_mag_phase(I1, I2, out)
% Task II: rebuild images from mixed magnitude and phase.
% Both images are the same size, so spectra can be swapped.
F1 = cfft(I1);  F2 = cfft(I2);
M1 = abs(F1);   P1 = angle(F1);
M2 = abs(F2);   P2 = angle(F2);

rand_phase  = 2*pi*rand(size(P1)) - pi;
noisy_phase = P1 + 0.8 * randn(size(P1));
flat_mag    = mean(M1(:)) * ones(size(M1));
k = fspecial('gaussian', 31, 8);
smooth_mag  = imfilter(M1, k, 'replicate');

zero_phase = fftshift(real(ifft2(ifftshift(M1))));

recs = { ...
    I1,                                   'Original image 1'; ...
    log(1 + abs(zero_phase)),             'Magnitude only (zero phase, log)'; ...
    icfft(M1 .* exp(1i * rand_phase)),    'Original magnitude + random phase'; ...
    icfft(M1 .* exp(1i * noisy_phase)),   'Original magnitude + noisy phase'; ...
    icfft(flat_mag .* exp(1i * P1)),      'Phase only (flat magnitude)'; ...
    icfft(smooth_mag .* exp(1i * P1)),    'Original phase + smoothed magnitude'; ...
    icfft(M2 .* exp(1i * P1)),            'Magnitude of img2 + phase of img1'; ...
    icfft(M1 .* exp(1i * P2)),            'Magnitude of img1 + phase of img2'};

fprintf('\nMagnitude/phase reconstruction, correlation with image 1 / image 2:\n');
for k = 3:size(recs, 1)
    R = recs{k, 1};
    fprintf('  %-38s %6.3f  %6.3f\n', recs{k, 2}, corr2(R, I1), corr2(R, I2));
end
show_grid(recs(:, 1)', recs(:, 2)', 2, 4, out, 'T2_mag_phase_reconstruction');
show_grid({I1, I2}, {'Image 1 (gray)', 'Image 2 (gray)'}, 1, 2, out, 'T2_mag_phase_inputs');
end
