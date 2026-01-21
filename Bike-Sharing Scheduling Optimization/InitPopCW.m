function Chrom = InitPopCW(NIND, N, cusnum, init_vc)
    Chrom = zeros(NIND, N);
    chrom = change(init_vc, N, cusnum);
    for j = 1:NIND
        Chrom(j,:) = chrom;
    end
end