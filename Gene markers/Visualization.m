function Visualization(cell_gene,n,gene_name,cell_name)
C=zeros(3,10*n);
C(1,:)=(1-cell_gene)*0.85;
C(2,:)=1-cell_gene;
C(3,:)=1;
for i=1:n
    for j=1:10
        X((i-1)*10+j,1)=j;
        Y((i-1)*10+j,1)=i;
    end  
end
clf
ZZZ=ones(4,10);
ZZZ=ZZZ+60;
image(ZZZ);
colormap(gray);
for i=1:10
line_shu(i,1)=i;
line_shu(i,2)=line_shu(i,1);
line_shu(i,3)=0.5;
line_shu(i,4)=n+0.5;
hold on;
plot(line_shu(i,1:2),line_shu(i,3:4),'-w','LineWidth',1.5);
end
for i=1:n
line_he(i,1)=0.5;
line_he(i,2)=10.5;
line_he(i,3)=i;
line_he(i,4)=i;
hold on;
plot(line_he(i,1:2),line_he(i,3:4),'-w','LineWidth',2);
end
scatter(X, Y,25, C','filled', 'SizeData', 175);
set(gca,'YTick',[1:1:4]);%设置要显示坐标刻度
set(gca, 'xticklabel', gene_name);
set(gca, 'yticklabel', cell_name);
xlim(gca,[0.5 10.5]);
% xtl=get(gca,'XTickLabel'); 
% xt=get(gca,'XTick');     
% xtextp=xt;                     
% ytextp(1:10)=n+0.7; 
% text(xtextp,ytextp,xtl,'HorizontalAlignment','right','rotation',40); 
set(gca,'XTickLabelRotation',46);
% xlabel('Engel');
title('Darmanis');
end