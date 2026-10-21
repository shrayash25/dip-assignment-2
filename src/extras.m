function extras(I, t, out)
% Beyond the tasks: convolution theorem and padding, Laplacian sharpening,
% high-frequency emphasis, homomorphic filtering.
[M, N] = size(I);
D = dist_grid(M, N);
F = cfft(I);

% 1. Convolution theorem: product of DFTs = circular convolution.
%    Without padding, the borders wrap around. With padding, result = linear convolution.
k = 15;
h = fspecial('gaussian', 2*k + 1, 5);
g_lin = conv2(I, h, 'full');
g_pad = real(ifft2(fft2(I, M + 2*k, N + 2*k) .* fft2(h, M + 2*k, N + 2*k)));
g_cir = real(ifft2(fft2(I) .* fft2(h, M, N)));
wrap = abs(g_cir - g_lin(1:M, 1:N));
fprintf('Convolution theorem: max error padded = %.2e, unpadded = %.2e\n', ...
    max(abs(g_pad(:) - g_lin(:))), max(wrap(:)));
show_grid({g_lin(1:M, 1:N), g_cir, wrap}, ...
    {'Spatial conv2', 'DFT product, no padding', 'Wraparound error'}, 1, 3, out, [t '_X1_conv_padding']);

% 2. Laplacian in the frequency domain: H = -4*pi^2*(u^2+v^2)
[v, u] = meshgrid(((1:N) - floor(N/2) - 1) / N, ((1:M) - floor(M/2) - 1) / M);
Hlap = -4 * pi^2 * (u.^2 + v.^2);
lap = icfft(F .* Hlap);
lap = lap / max(abs(lap(:)));
sharp = min(max(I - lap, 0), 1);

% 3. High-frequency emphasis + histogram equalization
Hhfe = 0.5 + 0.75 * (1 - lowpass('gauss', D, 40));
hfe = mat2gray(icfft(F .* Hhfe));
hfe_eq = histeq(hfe);

% 4. Homomorphic filtering: compress illumination, boost reflectance
gL = 0.5;  gH = 2;  c = 1;  D0 = 30;
Hh = (gH - gL) * (1 - exp(-c * D.^2 / D0^2)) + gL;
homo = mat2gray(exp(icfft(cfft(log(1 + I)) .* Hh)) - 1);

show_grid({I, lap, sharp, hfe, hfe_eq, homo}, ...
    {'Original', 'Laplacian (freq domain)', 'Laplacian sharpened', ...
     'High-frequency emphasis', 'HFE + histeq', 'Homomorphic'}, ...
    2, 3, out, [t '_X2_enhancement']);
end
