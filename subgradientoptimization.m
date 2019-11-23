clearvars
p11
%step 0% basic initializations
numnodes = dimX * dimY * 2;
pi = zeros(numnodes,1);
okcom = zeros(k,2);
directlink = zeros(k,1);
com = sortrows(com,1);
ergodicscom = com;
pov = 0;
bestpov = 0;
atiteration = 0;
h_history = [];

%step 1% Algorithm specific initializations. can be merged with step 0. 
lagrangianmult = zeros(numnodes,1);
ergodics = zeros(numnodes,1);
stepsizeoffset = 2;
lagrangianmult(1:numnodes) = 1/numnodes;

%step 2 Calculation for h_u and iterations begin here

for(iter = 1:1000)
    %Stepsize offset change after 10 iterations
    if(mod(iter,10) == 0)
        stepsizeoffset = 0.95*stepsizeoffset;
    end

%Initializing the values of pi = lagrangian mult for every iteration    
for(i=1:numnodes)
    pi(i) = lagrangianmult(i);
end

%GSP for initial path for independent k-pairs is returned in list
% if(iter == 1)
    list = gsp(dimX,dimY,pi,k,com);
    newlist = list;
% else
%     list = gsp(dimX,dimY,pi,numel(okcom)/2,okcom)
% end


% last = 0;
% for (i = 1 : k)
%     first = last+1;
%     slask = find(list(last+1:length(list)) == com(i,1));
%     last = slask(1)+first-1;
%     if (sum(pi(list(first:last))) < 1)
%         okcom = [okcom i]; 
%         newlist = [newlist; list(first:last)];
%     end
% end

%Algorithm for okcom -- okcom captures values of com when flag < 1 and sets
%others as -1 for loop uniformity. This may need changes as in few
%iterations all the com values maynot be present
okcom = com;
j = 1;
flag = 0;
okcompos = 1;
for(i = 1:k)
    while(com(i)~=list(j))
        flag = flag + pi(list(j));
        j = j+1;
    end
    flag = flag + pi(list(j));
    %fprintf('The value of Okcom for pair no %d is %f\n',i,flag);
    if(flag < 1)
        okcom(okcompos) = com(i);
        numokcom = numel(okcom);
        okcom(okcompos+numokcom/2) = com(i+k);
        okcompos = okcompos+1;
    else
        okcom(okcompos) = -1;
        numokcom = numel(okcom);
        okcom(okcompos+numokcom/2) = -1;
        okcompos = okcompos + 1;
    end
    j = j+1;
    flag = 0;
end

%Removal of -1's
numokcom = numel(okcom);
i = 1;
while(i<=(numokcom/2))
    if(okcom(i) == -1)
        okcom(i,:) = [];
        numokcom = numokcom-2;
        i = 1;
    else
        i = i+1;
    end        
end

%NewList creation from list s.t. for only okcom indecies. - working fine
j = 1;
i = 1;
loopcount = 1;
flag = 1;
flag2 = 1;
flag3 = 1;
numcom = numel(com);
for(loopcount = 1:numcom/2)
    flag = j; 
    while(com(i)~=list(j))
        j = j+1;
    end
    if(com(i) == okcom(flag2))
        j = flag;
        while(com(i)~=list(j))
            newlist(flag3) = list(j);
            j = j+1;
            flag3 = flag3+1;
        end
        newlist(flag3) = list(j);
        flag3 = flag3 + 1;
        flag2 = flag2 + 1;
    end
    j = j+1;
    i = i+1;
end

%Remoal of extra elements if any from prev iteration
if(flag3<=numel(newlist))
    newlist(flag3:numel(newlist)) = [];
end

%Calculation of h_u begins
%Hypothetical assignment of directlinks
flag = 0;
directlink = zeros(k,1);
for(i = 1:(numel(okcom)/2))
    flag = okcom(i);
    directlink(flag) = 1;
end
 
%calculation for h_mu
lagsum = 0;
for(i = 1:numnodes)
    lagsum = lagsum + lagrangianmult(i);
end

forh = zeros(k,1);
j=1;
for(i = 1:k)
    if(j<=numel(newlist))
        while(okcom(i)~=newlist(j))
            forh(i) = forh(i) + lagrangianmult(newlist(j));
            j = j+1;
        end
    end         
    %forh(i) = directlink(i) - forh(i);
    forh(i) = 1 - forh(i);
    j = j+1;
end
maxforh = sum(forh);
h = lagsum + maxforh;

%step 3%
%Adding links presence for all occurrences of node i.e outgoing and incoming
subgrad = zeros(numnodes,1);
for(i = 1:numnodes)
    flag = 0;
    for(j = 1:numel(newlist))
        if(i == newlist(j))
            flag = flag + 1;
        end
    end
    subgrad(i) = flag;
end

%Removing the outgoing occurences
for(i = 1:numnodes)
    for(j = 1:(numel(okcom))/2)
        if(i == okcom(j))
            subgrad(i) = subgrad(i) - 1;
        end
    end
end

% calculating subgrad i.e. 1-summation,summation for all i's
for(i = 1:numnodes)
    subgrad(i) = 1 - subgrad(i);
end

%calculating total subgradient for step length. subgrad squared and added
totsubgrad = 0;
for(i = 1:numnodes)
    totsubgrad = totsubgrad + subgrad(i) * subgrad(i);
end

%Calculating step length
steplen = stepsizeoffset*h/totsubgrad;
%fprintf('Value of steplen is %f\n',steplen);

%Calculating lagrangian multiplier
for(i=1:numnodes)
    lagrangianmult(i) = lagrangianmult(i) - steplen*subgrad(i);
    if(lagrangianmult(i) < 0)
        lagrangianmult(i) = 0;
    end
    if(lagrangianmult(i)~=0)
        %fprintf('value of lagrangian mult for %d is %f\n',i,lagrangianmult(i));
    end
end

h_history = [h_history, h];

%Step 4 - Ergodics
%Calculation of ergodics pi values
checkcount = zeros(dimX,1);
ergodicscom = com;
rule = 4;
summation1 = 0;
summation2 = 0;
if(iter>=2)
    for(i=1:numnodes)
        for(s = 0:iter-2)
            summation1 = summation1 + (s+1)^rule;
        end
        for(s = 0:iter-1)
            summation2 = summation2 + (s+1)^rule;
        end
        firstterm = (summation1/summation2)*ergodics(i);
        secondterm = ((iter^rule)/summation2)*lagrangianmult(i);
        ergodics(i) = firstterm + secondterm;
    end
    
%     Best lower bound - Heuristic approach
%     Using ergodics pi value, we find route for all pairs
    ergodicslist = gsp(dimX,dimY,ergodics,k,com);
    
%     Counting the num of nodes(outgoing nodes) each pair is taking up
    count = zeros(dimX,1);
    j = 1;
    blockposition = zeros(dimX,1);
    for(i=1:k)
        if(j <= numel(ergodicslist))
            blockposition(com(i)) = j;
            while(com(i) ~= ergodicslist(j))
                j = j+1;
                count(com(i)) = count(com(i)) + 1;
            end
        end
        j = j+1;
    end
%     Finding the Contact pair which takes maximum nodes in sequence - We should change
%     this to contact pair with maximum common pairs
    for(longiter = 1:k)
    for(i = 1:numel(find(count))>0)
        [maxcount, maxcountposition] = max(count);
        if(maxcount > 0)
        j = blockposition(maxcountposition);
        position = find(com == maxcountposition);
        while(com(position) ~= ergodicslist(j))
            for(check = 1:numel(ergodicslist))
                if(ergodicslist(j) == ergodicslist(check))
                    checkcount(maxcountposition) = checkcount(maxcountposition) + 1;
                end
            end
            j = j+1;
        end
        checkcount(maxcountposition) = checkcount(maxcountposition) - count(com(position));
        count(maxcountposition) = -1*count(maxcountposition);
        end
    end
%     Removal of contact pair if there are sharing nodes
    [maxcheckcount, maxcheckcountposition] = max(checkcount);
    if(maxcheckcount >= 1)
        j = blockposition(maxcheckcountposition);
        checkcount(maxcheckcountposition) = 0;
        count(maxcheckcountposition) = 0;
        position = find(com == maxcheckcountposition);
        while(com(position) ~= ergodicslist(j))
            ergodicslist(j) = -1;
            j = j+1;
        end
        ergodicslist(j) = -1;
        position = find(ergodicscom == maxcheckcountposition);
        ergodicscom(position , :) = [];
    end
    for(i = 1:numel(count))
        count(i) = -1*count(i);
        checkcount(i) = 0;
    end
    end
    ergodicslist(ergodicslist == -1) = [];
    
    
    % Obtain best LBD for each iteration
    % We have all the unique routes for this ergodic list of values. To
    % avoid these nodes completely we opt to choose these nodes values as
    % infinite
    %visagrid(dimX,dimY,ergodicslist,ergodicscom,ergodics,25);
    lbdpi = ergodics;
    for(i = 1:numel(ergodicslist))
        lbdpi(ergodicslist(i)) = numnodes;
    end
    for(i = 1:numel(com)/2)
        if(com(i) ~= ergodicscom(1:numel(ergodicscom)/2))
            reusablecom = [com(i),com(i+k)];
            reusablelist = gsp(dimX,dimY,lbdpi,1,reusablecom);
            if(numel(reusablelist)> 2 && sum(lbdpi(reusablelist)) < numnodes)
                ergodicscom = [ergodicscom; reusablecom];
                ergodicslist(numel(ergodicslist)+1:numel(ergodicslist)+numel(reusablelist)) = reusablelist(1:numel(reusablelist));
                for(j = 1:numel(reusablelist))
                    lbdpi(reusablelist(j)) = numnodes;
                end
            end
        end
    end   
    %visagrid(dimX,dimY,ergodicslist,ergodicscom,ergodics,25);
           
%     Finding the max number of contact pairs that fit in all iterations
    if(iter == 2)
        noc = numel(ergodicscom);
        maxconpov = numel(ergodicslist);
        maxconcom = ergodicscom;
        maxconlist = ergodicslist;
        maxconpi = ergodics;
        atiteration = iter;
    end
    if(noc <= numel(ergodicscom))
        if(noc == numel(ergodicscom))
            if(numel(ergodicslist) < maxconpov)
                maxconpov = numel(ergodicslist);
                maxconcom = ergodicscom;
                maxconlist = ergodicslist;
                maxconpi = ergodics;
                atiteration = iter;                
            end
        else
            noc = numel(ergodicscom);
            maxconpov = numel(ergodicslist);
            maxconcom = ergodicscom;
            maxconlist = ergodicslist;
            maxconpi = ergodics;
            atiteration = iter;
        end
    end
end
%End of iterations

fprintf('iter = %d\n',iter);
end
fprintf('Found at iteration number = %d\n',atiteration);
plot(h_history*4);
%Plotting to check
shift = 25;
% visagrid(dimX,dimY,newlist,okcom,pi,shift);
visagrid(dimX,dimY,maxconlist,maxconcom,maxconpi,shift);
