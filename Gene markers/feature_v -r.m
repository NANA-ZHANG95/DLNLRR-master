function gene_name=feature_v(Z,X,label,Genes,cell_name)
    n_space = length(unique(label));
    [num,level] = SIMLR_Feature_Ranking(Z,X');
    
    X_select=X(num(1:10),:);
    cell_level=zeros(10,n_space);
    for qq=1:n_space
        cell_level(:,qq)=cell_level(:,qq)+sum(X_select(:,find(label==qq)),2);
    end
    cell_level = mapminmax(cell_level,0,1);%归一化
    gene_name=Genes(num(1:10));   
    for i=1:n_space
        cell_level_heng((i-1)*10+1:i*10)=cell_level(1:10,i);
    end
    % 计算单细胞的稀疏度函数，大小代表稀疏性
    spare1=[];
    for i  = 1:length(cell_name)
        for j = 1:10
            X1 = X_select(j,find(label == i));
            spare1(10*(i-1)+j) = length(find(X1>0))/length(X1);
        end    
    end
    Visualization(cell_level_heng,n_space,gene_name,cell_name,spare1);
%%     figure 泡泡生成图
%     x1 = ones(1,5);
%     y1 = 1:5;
%     z1 = 0.2:0.2:1;
%     for i = 1:5
%         hold on
%         scatter(x1(i),y1(i),'k','filled', 'SizeData', 175*z1(i));
%     end
%     ylim([0,6])
%     xlim([0.5,1.5])

%%  Colorbar生成图
% A=1：10000；
% imagesc(A');
end