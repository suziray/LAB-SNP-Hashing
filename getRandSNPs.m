function SNPs = getRandSNPs(N, n) %N samples, n SNPs

load('distribution');
[~,s] = size(chr005);

fs = round(rand(1, n) * s, 0);
for i = 1:n
    fs(i) = chr005(fs(i));
end

SNPs = zeros(N, 2*n);
for i = 1:n
    SNPs(:, 2*i-1:2*i) = randSNPGenerator(N,fs(i));
end