function I = icfft(F)
% Inverse of cfft, real part only.
I = real(ifft2(ifftshift(F)));
end
