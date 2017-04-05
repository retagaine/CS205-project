fid=fopen('BFS_2.txt','r');
n = 122;
P = zeros(n,5); %1029258
fgetl(fid);
for ii = 1:n
    for jj = 1:5
        P(ii,jj) = fscanf(fid,'%f',1);
    end
    fgetl(fid);
end
fclose(fid);
[Q,I] = sort(P(:,4)); 
B = P(I,:);
col = colormap;
%hold on;
nd = 0.308473;
a = 1.0104;
c = 1.0074;
spos_Si = [0.33333333,  0.66666667,  0.93750000-1; 0.00000, 0.00000, 0.1875; 0.66666667,  0.33333333,  0.43750000; 0.000000,  0.000000,  0.68750000;  0.33333333,  0.66666667,  0.93750000; 0.00000, 0.00000, 0.1875+1];
cell2 = [nd, 0*a, 0*a; -nd/2,nd/2*sqrt(3), 0*a; 0*1.0, 0*1.0, 10.086*c*nd/3.078*a/2.57218587467527*2.51866888630220];
yval = (spos_Si(3,:)-spos_Si(2,:))*cell2;
yval = 2*yval(2);
%h = histogram2(P(:,1),P(:,2),[5,5],'DisplayStyle','tile','ShowEmptyBins','on');
t = -1;
maxer = B(1,5);
for jj = 1:n
    hold on;
    if B(jj,4) ~= t
        xlim([-0.8,0.8]);
        ylim([-0.8,0.8]);
        axis square
        colorbar('TickLabels',round(100*linspace(0,maxer,11)/10^round(log(maxer)/log(10)))/100*10^round(log(maxer)/log(10)));
        set(gcf,'Color','w');
        xlabel('\Delta x (nm)');
        ylabel('\Delta y (nm)');
        title(['Probability map at: ',num2str(t,'%.4f'),' s (1300 K)'])
        t = B(jj,4);
        box on;
        figure
        hold on;
        maxer = max(B(jj:end,5));
        for ii = 1:jj
            scatter(B(ii,1),B(ii,2),'s','filled','MarkerEdgeColor',col(round(63*min([1,B(ii,5)/maxer]))+1,:),'MarkerFaceColor',col(round(63*min([1,B(ii,5)/maxer]))+1,:),'SizeData',800)
        end
        for ii = jj+1:n
           scatter(B(ii,1),B(ii,2),'s','filled','MarkerEdgeColor',col(round(63*min([0,B(ii,5)/maxer]))+1,:),'MarkerFaceColor',col(round(63*min([0,B(ii,5)/maxer]))+1,:),'SizeData',800)
        end
    end
    scatter(B(jj,1),B(jj,2),'s','filled','MarkerEdgeColor',col(round(63*min([1,B(jj,5)/maxer]))+1,:),'MarkerFaceColor',col(round(63*min([1,B(jj,5)/maxer]))+1,:),'SizeData',800)
end
xlim([-0.8,0.8]);
ylim([-0.8,0.8]);
axis square
colorbar('TickLabels',round(100*linspace(0,maxer,11)/10^round(log(maxer)/log(10)))/100*10^round(log(maxer)/log(10)));
set(gcf,'Color','w');
xlabel('\Delta x (nm)');
ylabel('\Delta y (nm)');
box on;
title(['Probability map at: ',num2str(t,'%.4f'),' s (1300 K)'])  