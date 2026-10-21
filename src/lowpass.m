function H = lowpass(type, D, D0, n)
% Low-pass transfer function. type: 'ideal', 'gauss', 'butter' (order n).
if nargin < 4, n = 2; end
switch type
    case 'ideal',  H = double(D <= D0);
    case 'gauss',  H = exp(-(D.^2) ./ (2 * D0^2));
    case 'butter', H = 1 ./ (1 + (D ./ D0).^(2 * n));
end
end
