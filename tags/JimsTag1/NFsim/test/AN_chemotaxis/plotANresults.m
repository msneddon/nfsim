
%%%% Add the path to the NFanaylsis tools
addpath '../../tools/NFanalysis/';

%%%% Plot the average methylation level
x=tblread('an2_nf.gdat','\t');
meanMethLevel=(0.*x(:,2)+ 1.*x(:,3)+ 2.*x(:,4)+3.*x(:,5)+ 4.*x(:,6)+ 5.*x(:,7)+ 6.*x(:,8)+ 7.*x(:,9)+ 8.*x(:,10))./sum(x(:,2:10),2);
figure; plot(x(:,1)./60,meanMethLevel,'LineWidth',3,'Color','b'); hold on; 
axis([0,x(end,1)./60,0,8]);
xlabel('Time (minutes)','FontName', 'Arial','fontSize',18);
ylabel('Average Methylation Level','FontName', 'Arial','fontSize',18);

% make things nice
grid on;
set(gca, 'TickDir', 'out');
set(gcf, 'color', 'white');
set(gca, 'FontName', 'Arial');
set(gca, 'fontSize', 15);


%%%%% Plot a surface of how the methylation distribution changes in time
figure;
surf([0:8],x(:,1)./60,x(:,2:10)./sum(x(1,2:10)),'EdgeColor','none','FaceAlpha',1,'FaceLighting','phong');
hold on;
grid on;
ylabel('Time (minutes)','FontName', 'Arial','fontSize',18);
xlabel('Methylation Level','FontName', 'Arial','fontSize',18);
zlabel('Fraction of Receptors','FontName', 'Arial','fontSize',18);
set(gca, 'TickDir', 'out');
set(gcf, 'color', 'white');
set(gca, 'FontName', 'Arial');
set(gca, 'fontSize', 15);



%%%% Plot the average activity of each receptor dimer
r=readNFdump();
pOn=getAvgLocalFunctionValue(r,'RD');
figure; plot(getTimeArray(r)./60,pOn,'LineWidth',3,'Color','b'); hold on; 
axis([0,max(getTimeArray(r)./60),0,0.5]);
xlabel('Time (minutes)','FontName', 'Arial','fontSize',18);
ylabel('Average Methylation Level','FontName', 'Arial','fontSize',18);

grid on;
set(gca, 'TickDir', 'out');
set(gcf, 'color', 'white');
set(gca, 'FontName', 'Arial');
set(gca, 'fontSize', 15);


%%%% Remove the 
%rmpath '../../tools/NFanalysis/';