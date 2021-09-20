%%This script visualize the scale of "individual" eigenfunction
figure;
EF_no=2;%This is the eigenfunction numer for visualization. You can change to another number
patch(leftSurface,'FaceVertexCData',eigenFunctionsL(:,EF_no),'EdgeColor','none','FaceColor','interp');
axis off;
axis image;
hold on 
patch(rightSurface,'FaceVertexCData',eigenFunctionsR(:,EF_no),'EdgeColor','none','FaceColor','interp');

view([0,90]); %view([270,0]) for left himisphere only; view([90 0]) for RH only
set(gcf,'color','w');
