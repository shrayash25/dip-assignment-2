function H = bandreject(type, D, C0, W, n)
% Band-reject transfer function, center radius C0, width W.
% Band-pass is 1 - bandreject(...).
if nargin < 5, n = 2; end
switch type
    case 'ideal',  H = double(~(D >= C0 - W/2 & D <= C0 + W/2));
    case 'gauss',  H = 1 - exp(-((D.^2 - C0^2) ./ (D .* W + eps)).^2);
    case 'butter', H = 1 ./ (1 + ((D .* W) ./ (D.^2 - C0^2 + eps)).^(2 * n));
end
end
