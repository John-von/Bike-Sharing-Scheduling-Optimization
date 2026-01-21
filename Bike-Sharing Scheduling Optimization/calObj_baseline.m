function ObjV = calObj_baseline(Chrom, cusnum, cap, demands, a, b, L, s, dist, alpha, beta)
    NIND = size(Chrom, 1);
    ObjV = zeros(NIND, 1);
    for i = 1:NIND
        [VC, ~, ~, ~, ~] = decode_baseline(Chrom(i,:), cusnum, cap, demands, a, b, L, s, dist);
        ObjV(i) = costFuction(VC, a, b, s, L, dist, demands, cap, alpha, beta);
    end
end