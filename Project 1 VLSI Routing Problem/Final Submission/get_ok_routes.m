function [ ok_routes, new_nodelist ] = get_ok_routes(pi, k, routes, nodelist)
    last = 0;
    ok_routes = [];
    new_nodelist = [];
    for i = 1:k
        first = last + 1;
        slask = find(nodelist(last+1:length(nodelist)) == routes(i, 1));
        last = slask(1) + first - 1;
        if (sum(pi(nodelist(first:last))) < 1)
            ok_routes = [ok_routes; nodelist(last)];
            new_nodelist = [new_nodelist; nodelist(first:last)];
        end
    end
end
