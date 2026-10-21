function I = load_gray(file, n)
% Read image as gray double in [0,1], center-crop to square, resize to n x n.
I = imread(file);
if size(I, 3) == 3, I = rgb2gray(I); end
I = im2double(I);
[h, w] = size(I);
s = min(h, w);
r0 = floor((h - s) / 2);  c0 = floor((w - s) / 2);
I = I(r0+1:r0+s, c0+1:c0+s);
I = imresize(I, [n n]);
I = min(max(I, 0), 1);
end
