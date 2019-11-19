clearvars
p6
%step 0%
numnodes = dimX * dimY * 2;
pi = zeros(numnodes,1);
okcom = zeros(k,2);
directlink = zeros(k,1);
%newlist = zeros(numnodes,1);

%step 1%
lagrangianmult = zeros(numnodes,1);
stepsizeoffset = 2;

%step 2 for one iteration

for(iter = 1:50)
    if(mod(iter,10) == 0)
        stepsizeoffset = 0.95*stepsizeoffset;
    end

for(i=1:numnodes)
    pi(i) = lagrangianmult(i);
end

if(iter == 1)
    list = gsp(dimX,dimY,pi,k,com)
    newlist = list;
else
    list = gsp(dimX,dimY,pi,k,okcom)
end
% listcount = numel(list)
% count = zeros(k,1)
% j = 1;
% for(i=1:k)
%     if(j<=numel(list))
%         while(i~=list(j))
%             j = j+1;
%             count(i) = count(i) + 1;
%         end
%     end
%     count(i) = count(i) - 1;
%     fprintf('count of %d = %d\n',i,count(i));
%     j = j+1;
% end
% for(i = 1:k)
%     if(count(i)*0.1>=1)
%         count(i) = 0;
%     end
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

%for okcom and newlist
j = 0;
flag = 0;
okcompos = 1;
for(i = 1:k)
    j = j+1;
    while(i~=list(j))
        flag = flag + pi(list(j));
        j = j+1;
    end
    fprintf('The value for okcom is %f\n',flag);
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
end
numokcom = numel(okcom);
i = 1;
while(i<=(numokcom/2))
    if(okcom(i) == -1)
        okcom(i,:) = [];
        numokcom = numokcom-2;
        i = i-1;
    else
        i = i+1;
    end        
end

j = 1;
flag = 1;
flag2 = 1;
flag3 = 1;
for(i = 1:k)
    flag = j;
    while(i~=list(j))
        j = j+1;
    end
    if(i == okcom(flag2))
        j = flag;
        while(i~=list(j))
            newlist(flag3) = list(j);
            j = j+1;
            flag3 = flag3+1;
        end
        newlist(flag3) = list(j);
        flag3 = flag3 + 1;
        flag2 = flag2 + 1;
    end
    j = j+1;
end


flag = 0;
for(i = 1:(numel(okcom)/2))
    flag = okcom(i)
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
    forh(i) = directlink(i) - forh(i)
    j = j+1;
end
maxforh = max(forh);
h = lagsum + maxforh;

%step 3%
%Adding for all occurrences of node i.e outgoing and incoming
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
    subgrad(i) = 1 - subgrad(i)
end
%calculating total subgradient for step length. subgrad squared and added
totsubgrad = 0;
for(i = 1:numnodes)
    totsubgrad = totsubgrad + subgrad(i) * subgrad(i);
end
%Calculating step length
steplen = stepsizeoffset*h/totsubgrad;
fprintf('Value of steplen is %f\n',steplen);

%Calculating lagrangian multiplier
for(i=1:numnodes)
    lagrangianmult(i) = lagrangianmult(i) - steplen*subgrad(i);
    if(lagrangianmult(i) < 0)
        lagrangianmult(i) = 0;
    end
    if(lagrangianmult(i)~=0)
        fprintf('value of lagrangian mult is %f\n',lagrangianmult(i));
    end
end
%End of iterations
end
%Plotting to check
shift = 25;
visagrid(dimX,dimY,newlist,okcom,pi,shift);