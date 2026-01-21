function SelCh = LocalSearch_hybrid(SelCh, cusnum, cap, demands, a, b, L, s, dist, alpha, beta, operator_idx)
    D = 15;
    toRemove = 15;
    [row, N] = size(SelCh);
    for i = 1:row
        [VC, ~, ~, ~, ~] = decode_hybrid(SelCh(i,:), cusnum, cap, demands, a, b, L, s, dist, operator_idx);
        CF = costFuction(VC, a, b, s, L, dist, demands, cap, alpha, beta);
        [removed, rfvc] = Remove(cusnum, toRemove, D, dist, VC);
        [ReIfvc, ~] = Re_inserting(removed, rfvc, L, a, b, s, dist, demands, cap);
        RCF = costFuction(ReIfvc, a, b, s, L, dist, demands, cap, alpha, beta);
        if RCF < CF
            chrom = change(ReIfvc, N, cusnum);
            if length(chrom) == N
                SelCh(i,:) = chrom;
            end
        end
    end
end