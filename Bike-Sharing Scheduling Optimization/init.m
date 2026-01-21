function [init_vc] = init(cusnum, a, demands, cap)
    j = ceil(rand * cusnum);
    k = 1;
    init_vc = cell(k, 1);
    if j == 1
        seq = 1:cusnum;
    elseif j == cusnum
        seq = [cusnum, 1:j-1];
    else
        seq1 = 1:j-1;
        seq2 = j:cusnum;
        seq = [seq2, seq1];
    end
    route = [];
    load = 0;
    i = 1;
    while i <= cusnum
        if load + demands(seq(i)) <= cap
            load = load + demands(seq(i));
            if isempty(route)
                route = [seq(i)];
            elseif length(route) == 1
                if a(seq(i)) <= a(route(1))
                    route = [seq(i), route];
                else
                    route = [route, seq(i)];
                end
            else
                lr = length(route);
                flag = 0;
                for m = 1:lr-1
                    if (a(seq(i)) >= a(route(m))) && (a(seq(i)) <= a(route(m+1)))
                        route = [route(1:m), seq(i), route(m+1:end)];
                        flag = 1;
                        break
                    end
                end
                if flag == 0
                    route = [route, seq(i)];
                end
            end
            if i == cusnum
                init_vc{k,1} = route;
                break
            end
            i = i + 1;
        else
            init_vc{k,1} = route;
            route = [];
            load = 0;
            k = k + 1;
        end
    end
end