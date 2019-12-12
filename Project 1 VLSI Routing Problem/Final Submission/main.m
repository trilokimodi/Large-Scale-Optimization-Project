tic
p10
num_nodes = dimX * dimY * 2;
pi = ones(dimX*dimY*2, 1) / num_nodes;
theta = 1.99999;
h_history = [];
shift = 25;
LBD = 1 / num_nodes;

best_LBD = -1;
best_okcom = [];
best_nodelist = [];
LBD_history = [];
sum_pi = zeros(dimX*dimY*2,1);
flag = 0;


for iteration = 1:1000
    % stepsize offset change after 10 iterations
    if (mod(iteration, 10) == 0)
        theta = 0.95 * theta;
        fprintf("iteration: %4d   theta: %1.3f\n", iteration, theta);
    end
    
    % initializing the values of pi = lagrangian multiplier for every iteration
    nodelist = gsp(dimX, dimY, pi, k, com);
    
    [okcom, new_nodelist] = get_ok_routes(pi, k, com, nodelist);
    [hokcom, hnew_nodelist] = heuristics(dimX, dimY, pi, okcom, new_nodelist);
    
    LBD = numel(hokcom);
    if (LBD > best_LBD)
        iter = iteration;
        best_LBD = LBD;
        best_okcom = hokcom;
        best_nodelist = hnew_nodelist;
    elseif(LBD == best_LBD)
        if(numel(best_nodelist) > numel(hnew_nodelist))
            iter = iteration;
            flag = flag + 1;
            best_LBD = LBD;
            best_okcom = hokcom;
            best_nodelist = hnew_nodelist;
        end
    end
    
    % calculation of h
    subh = zeros(numel(okcom), 1);
    sub_grad = zeros(num_nodes, 1);
    
    j = 1;
    for i = 1:numel(okcom)
        while (okcom(i) ~= new_nodelist(j))
            % add 1 count for all nodes used
            sub_grad(new_nodelist(j)) = sub_grad(new_nodelist(j)) + 1;
            
            % add the cost for all links used
            subh(i) = subh(i) + pi(new_nodelist(j));
            j = j + 1;
        end
        % add 1 count for the exit node
        sub_grad(new_nodelist(j)) = sub_grad(new_nodelist(j)) + 1;
        
        % add cost for the link used to the exit node
        subh(i) = subh(i) + pi(new_nodelist(j));
        
        subh(i) = 1 - subh(i);
        j = j + 1;
    end
    
    h = sum(pi) + sum(subh);
    
    % calculating subgrad i.e. 1 - summation, summation for all i's
    for i = 1:num_nodes
        sub_grad(i) = 1 - sub_grad(i);
    end
    
    % calculating total subgradient for step length. subgrad squared and added
    tot_subgrad = sum(sub_grad.*sub_grad);
    
    % calculating step length
    step_length = (theta * (h - best_LBD)) / tot_subgrad;
    
    % calculating lagrangian multiplier
    for i = 1:num_nodes
        pi(i) = max(0, pi(i)-step_length*sub_grad(i));
        sum_pi(i) = sum_pi(i) + pi(i);
    end
    
    h_history = [h_history, h];
    LBD_history = [LBD_history; best_LBD];
end
toc
visagrid(dimX, dimY, best_nodelist, best_okcom, pi, shift);
toc