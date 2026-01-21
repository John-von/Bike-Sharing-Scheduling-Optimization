function dChrom = deal_Repeat(Chrom)
    N = size(Chrom, 1);
    len = size(Chrom, 2);
    dChrom = unique(Chrom, 'rows');
    Nd = size(dChrom, 1);
    newChrom = zeros(N-Nd, len);
    for i = 1:N-Nd
        newChrom(i,:) = randperm(len);
    end
    dChrom = [dChrom; newChrom];
end