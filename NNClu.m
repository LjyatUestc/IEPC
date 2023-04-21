function clustering = NNClu(data,cnum,redataind,seeds)
redataind(seeds) = [];
sim = pdist2(data,data,'euclidean');
flag = 1;
toc;
clustering = cell(cnum,1);
for i = 1:cnum
    clustering{i,1} = seeds(i);%%%初始化聚类种子点
end
while ~isempty(redataind)
    if flag == 1
    rer = size(redataind,1);
    minval = [];
    label = [];    
    for rd = 1:rer
        td = [];
        clulabel = [];
        for i = 1:cnum      
            C = clustering{i,1};
            rc = size(C,2);
            for j = 1:rc
                d = sim(redataind(rd),C(j));
                td = [td d];
                clulabel = [clulabel i];
            end     
        end
        [val ind] = min(td);
        minval = [minval val];
        label = [label clulabel(ind)];
    end
    [val ind] = min(minval);
    explabel = label(ind);
    Ctemp=clustering{explabel,1};
    Ctemp(end+1)=redataind(ind);
    clustering{explabel,1}=Ctemp;
    temp = redataind(ind);
    redataind(ind) = [];
    minval(ind) = [];
    label(ind) = [];
    flag = 0;
    else
        rer = size(redataind,1);
        for m = 1:rer
            td = sim(redataind(m),temp);
            if td < minval(m)
                minval(m) = td;
                label(m) = explabel;
            end             
        end
        [val ind] = min(minval);
        explabel = label(ind);
        Ctemp=clustering{explabel,1};
        Ctemp(end+1)=redataind(ind);
        clustering{explabel,1}=Ctemp;
        temp = redataind(ind);
        redataind(ind) = [];
        minval(ind) = [];
        label(ind) = [];
    end
end