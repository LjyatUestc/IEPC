function [ind,dy,rho,dypro] = IEPks(Z,knn)
parpool(4);
Md1 = createns(Z);
dy = [];
rows = size(Z,1);
knnindex = zeros(rows,knn-1);

parfor i = 1 : rows  %%%%%%%%%%%%构建kd树求k最近邻
    [midx,md] = knnsearch(Md1,Z(i,:),'K',knn);
    knnindex(i,:) = midx(2:end);  
end
parfor j = 1 : knn-1 %%%%%%%%%%%%%%%%%%统计动态增长逆k最近邻
    det = tabulate(knnindex(:,j));
    incre = [det(:,2);zeros(rows-length(det(:,2)),1)];
%     index = det(:,1);
    dy = [dy incre];
%     dy(index,j) = det(:,2);
end
delete(gcp('nocreate'));
dyentro = dy;
dypro = dy;
for j = 1 : knn-1
    te = tabulate(dy(:,j));
    pro = te(:,3)./100;
    entro = -log2(pro).*(pro);
    leng = size(te(:,1),1);
    for jj = 1:leng
        in = find(dyentro(:,j) == te(jj,1));
        dyentro(in,j) = entro(jj);
        dypro(in,j) = pro(jj);
    end
    
end
% dyentro(:,knn) = 1./(sum(dyentro,2)./sum(dy,2));%%%%计算动态信息熵
dyentro(:,knn) = sum(dyentro,2).*sum(dy,2);
[rho, ind] = sort(dyentro(:,end),'descend');
% ttt = 0;
