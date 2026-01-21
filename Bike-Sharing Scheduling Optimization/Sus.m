function NewChrIx = Sus(FitnV, Nsel)
    [Nind, ~] = size(FitnV);
    if Nind == 0 || Nsel <= 0
        NewChrIx = [];
        return;
    end
    if all(FitnV == 0) || all(isnan(FitnV)) || all(isinf(FitnV))
        NewChrIx = randperm(Nind, min(Nsel, Nind));
        return;
    end
    cumfit = cumsum(FitnV);
    if cumfit(end) == 0 || isnan(cumfit(end)) || isinf(cumfit(end))
        NewChrIx = randperm(Nind, min(Nsel, Nind));
        return;
    end
    trials = cumfit(Nind) / Nsel * (rand + (0:Nsel-1)');
    Mf = cumfit(:, ones(1, Nsel));
    Mt = trials(:, ones(1, Nind))';
    [NewChrIx, ~] = find(Mt < Mf & [zeros(1, Nsel); Mf(1:Nind-1, :)] <= Mt);
    if isempty(NewChrIx)
        NewChrIx = randperm(Nind, min(Nsel, Nind));
    elseif length(NewChrIx) < Nsel
        remaining = setdiff(1:Nind, NewChrIx);
        need = Nsel - length(NewChrIx);
        if ~isempty(remaining)
            extra = remaining(randperm(length(remaining), min(need, length(remaining))));
            NewChrIx = [NewChrIx; extra(:)];
        end
    elseif length(NewChrIx) > Nsel
        NewChrIx = NewChrIx(randperm(length(NewChrIx), Nsel));
    end
    if ~isempty(NewChrIx)
        [~, shuf] = sort(rand(length(NewChrIx), 1));
        NewChrIx = NewChrIx(shuf);
    end
end