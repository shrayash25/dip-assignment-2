function task2_synthetic(out)
% Task II: synthetic image with three frequency patterns, select with BPF / BRF.
N = 256;
[x, y] = meshgrid(0:N-1, 0:N-1);
f = [8 30 50*sqrt(2)];                           % radial frequency of each pattern
p1 = cos(2*pi*f(1)*x/N);                         % low:  vertical stripes
p2 = cos(2*pi*f(2)*y/N);                         % mid:  horizontal stripes
p3 = cos(2*pi*50*(x + y)/N);                     % high: diagonal stripes
S = (p1 + p2 + p3) / 3;
pats = {p1, p2, p3};

F = cfft(S);
D = dist_grid(N, N);
W = 10;  n = 4;

show_grid({p1, p2, p3, S, logmag(F)}, ...
    {sprintf('Pattern 1 (f=%.0f)', f(1)), sprintf('Pattern 2 (f=%.0f)', f(2)), ...
     sprintf('Pattern 3 (f=%.0f)', f(3)), 'Synthetic image', 'Spectrum'}, ...
    2, 3, out, 'T2_synthetic_input');

ims = {};  tt = {};
for k = 1:3
    Hbp = 1 - bandreject('butter', D, f(k), W, n);
    G = icfft(F .* Hbp);
    ims{end+1} = G;
    tt{end+1} = sprintf('BPF keep f=%.0f, MSE %.4f', f(k), mean((G(:) - pats{k}(:)/3).^2));
end
for k = 1:3
    Hbr = bandreject('butter', D, f(k), W, n);
    G = icfft(F .* Hbr);
    ref = S - pats{k}/3;
    ims{end+1} = G;
    tt{end+1} = sprintf('BRF reject f=%.0f, MSE %.4f', f(k), mean((G(:) - ref(:)).^2));
end
for k = 1:3
    ims{end+1} = 1 - bandreject('butter', D, f(k), W, n);
    tt{end+1} = sprintf('BPF mask f=%.0f', f(k));
end
show_grid(ims, tt, 3, 3, out, 'T2_synthetic_bpf_brf');
show_spectra([{S}, ims(1:3)], [{'Synthetic image'}, tt(1:3)], out, 'T2_synthetic_bpf_spectra');
show_spectra([{S}, ims(4:6)], [{'Synthetic image'}, tt(4:6)], out, 'T2_synthetic_brf_spectra');
end
