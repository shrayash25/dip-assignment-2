function p = psnr_db(X, ref)
% PSNR in dB for images in [0,1].
p = 10 * log10(1 / mean((X(:) - ref(:)).^2));
end
