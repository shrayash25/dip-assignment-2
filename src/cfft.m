function F = cfft(I)
% Centered 2-D DFT (DC in the middle).
F = fftshift(fft2(I));
end
