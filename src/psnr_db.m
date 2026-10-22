function p = psnr_db(X, ref)
% PSNR in dB for images in [0,1].
mse = mean((X(:) - ref(:)).^2);
if mse < 1e-20, p = Inf; else, p = 10 * log10(1 / mse); end
end
