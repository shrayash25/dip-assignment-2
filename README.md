# EEE F435 DIP – Assignment 2

## Run
From this folder:

```
matlab -batch "run_all"
```

Or open MATLAB in this folder and type `run_all`.

Needs MATLAB with the Image Processing Toolbox.

## Output
- `results/*.png` – every figure (original and processed images, spectra, plots)
- `results/log.txt` – printed numbers (PSNR, edge ratio, SNR, errors)

## Files
| File | Content |
|---|---|
| `run_all.m` | Entry point. Sets BITS ID (A = 3, B = 3), loads images, runs all tasks |
| `src/task1_properties.m` | Task I: spectrum, fftshift, shift, scale, rotation |
| `src/task2_mag_phase.m` | Magnitude / phase reconstruction |
| `src/task2_filters.m` | Ideal, Gaussian, Butterworth LPF/HPF, cutoffs, orders, noise |
| `src/task2_bands.m` | Low / mid / high band split |
| `src/task2_synthetic.m` | Synthetic 3-pattern image, BPF and BRF |
| `src/task2_periodic_noise.m` | Periodic noise, notch and band-reject |
| `src/task2_compression.m` | Spectrum truncation vs reconstruction error |
| `src/extras.m` | Extra: convolution theorem and padding, Laplacian, high-frequency emphasis, homomorphic |
| `images/` | Image 1 and Image 2 |
