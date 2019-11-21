%%
clc
clearvars
p10
com = sortrows(com,1);
num_nodes = dimX * dimY * 2;
pi = zeros(dimX * dimY * 2, 1);
lagrangian_multiplier = ones(dimX * dimY * 2, 1);
lagrangian_multiplier = lagrangian_multiplier * (1/num_nodes);
stepsizeoffset = 1.99999;
h_history = [];

for iteration = 1:1000
    % stepsize offset change after 10 iterations
    if(mod(iteration, 10) == 0)
        stepsizeoffset = 0.95*stepsizeoffset;
    end
    
    % initializing the values of pi = lagrangian multiplier for every iteration
    pi = lagrangian_multiplier;
    
    nl = gsp(dimX, dimY, pi, k, com);
    
    last = 0;
    okcom = [];
    newnl = [];
    for i = 1:k
        first = last+1;
        slask = find(nl(last+1:length(nl)) == com(i,1));
        last = slask(1)+first-1;
        if (sum(pi(nl(first:last))) < 1)
            okcom = [okcom i]; %#ok<AGROW>
            newnl = [newnl; nl(first:last)]; %#ok<AGROW>
        end
    end
    
    % calculation of h
    direct_link = zeros(k,1);
    direct_link(okcom) = 1;
    
    % TODO: use com and nl instead of okcom and newnl?
    % if not, then we don't need direct_link at all
    forh = zeros(k,1);
    j=1;
    for i = 1:k
        if(j<=numel(newnl))
            while(j<=numel(newnl) && okcom(i)~=newnl(j))
                forh(i) = forh(i) + lagrangian_multiplier(newnl(j));
                j = j+1;
            end
        end
        forh(i) = direct_link(i) - forh(i);
        j = j+1;
    end
    maxforh = max(forh);
    h = sum(lagrangian_multiplier) + maxforh;
    
    % step 3%
    % adding for all occurrences of node i.e outgoing and incoming
    subgrad = zeros(num_nodes,1);
    for i = 1:num_nodes
        for j = 1:numel(newnl)
            if(i == newnl(j))
                subgrad(i) = subgrad(i) + 1;
            end
        end
    end
    
    % removing the outgoing occurences
    for i = 1:num_nodes
        for j = 1:(numel(okcom))
            if(i == okcom(j))
                subgrad(i) = subgrad(i) - 1;
            end
        end
    end
    
    % calculating subgrad i.e. 1-summation,summation for all i's
    for i = 1:num_nodes
        subgrad(i) = 1 - subgrad(i);
    end
    
    % calculating total subgradient for step length. subgrad squared and added
    totsubgrad = 0;
    for i = 1:num_nodes
        totsubgrad = totsubgrad + subgrad(i) * subgrad(i);
    end
    
    % calculating step length
    step_length = (stepsizeoffset*h) / totsubgrad;
    
    % calculating lagrangian multiplier
    for i=1:num_nodes
        lagrangian_multiplier(i) = max(0, lagrangian_multiplier(i) - step_length*subgrad(i));
    end
    
    h_history = [h_history h]; %#ok<AGROW>
end

shift = 25;
visagrid(dimX,dimY,nl,com,pi,shift);

% include a plot of h values in the bottom
plot(h_history*4)
