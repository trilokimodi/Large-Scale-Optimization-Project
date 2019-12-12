function [ okcom, new_nodelist ] = heuristics( dimX, dimY, pi, okcom, ...
                                               new_nodelist)
% calculate nodes used
nodes_used = zeros(size(pi));
for i = 1:length(nodes_used)
    nodes_used(i) = length(find(new_nodelist == i));
end

% we can be lucky here and get a feasible solution without any heuristics
if (max(nodes_used) > 1)
    
    % calculate infeasible nodes
    infeasible_nodes = find(nodes_used >= 2);
    
    while (~isempty(infeasible_nodes))
        % find the "biggest offendertry" route
        % that is: the route with the most number of infeasible nodes
        max_infeasible_nodes_used = 0;
        last = 0; bad_last = 0; bad_first = 0;
        for i = 1:size(okcom)
            first = last+1;
            slask = find(new_nodelist(last+1:length(new_nodelist)) == okcom(i));
            last = slask(1)+first-1;
            
            infeasible_nodes_used = sum(ismember(new_nodelist(first:last), ...
                                                 infeasible_nodes));
            if infeasible_nodes_used > max_infeasible_nodes_used
                bad_first = new_nodelist(first);
                bad_first_index = first;
                
                bad_last = new_nodelist(last);
                bad_last_index = last;
            end
        end
        
        % remove biggest offender
        new_nodelist(bad_first_index:bad_last_index) = [];
        okcom(okcom == bad_last) = [];
        
        % try to find a feasible route
        % we set all used nodes to be expensive so that gsp avoids them
        expensive_pi = pi;
        expensive_pi(new_nodelist) = 1000;
        maybe_nl = gsp(dimX, dimY, expensive_pi, 1, [bad_last bad_first]);
        
        % if we found a feasible route, use it
        if length(new_nodelist) >= 2 && sum(expensive_pi(maybe_nl)) < 1
            okcom = [okcom; bad_last]; %#ok<AGROW>
            new_nodelist = [new_nodelist; maybe_nl]; %#ok<AGROW>
        end
        
        % update nodes used
        for i = 1:length(nodes_used)
            nodes_used(i) = length(find(new_nodelist == i));
        end
        
        % update infeasible nodes
        infeasible_nodes = find(nodes_used >= 2);
    end
end
end
