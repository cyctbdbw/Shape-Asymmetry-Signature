EF_num; %this is the total number of eigenfunctions/scales
coefsx_l=decomposeEigenFunction(eigenFunctionsL(:,1:n),leftSurface.vertices(:,1));
coefsy_l=decomposeEigenFunction(eigenFunctionsL(:,1:n),leftSurface.vertices(:,2));
coefsz_l=decomposeEigenFunction(eigenFunctionsL(:,1:n),leftSurface.vertices(:,3));

new_leftSurface.vertices(:,1) = eigenFunctionsL(:,1:n)*coefsx_l;
new_leftSurface.vertices(:,2) = eigenFunctionsL(:,1:n)*coefsy_l;
new_leftSurface.vertices(:,3) = eigenFunctionsL(:,1:n)*coefsz_l; 

coefsx_r=decomposeEigenFunction(eigenFunctionsR(:,1:n),rightSurface.vertices(:,1));
coefsy_r=decomposeEigenFunction(eigenFunctionsR(:,1:n),rightSurface.vertices(:,2));
coefsz_r=decomposeEigenFunction(eigenFunctionsR(:,1:n),rightSurface.vertices(:,3));

new_rightSurface.vertices(:,1) = eigenFunctionsR(:,1:n)*coefsx_r;
new_rightSurface.vertices(:,2) = eigenFunctionsR(:,1:n)*coefsy_r;
new_rightSurface.vertices(:,3) = eigenFunctionsR(:,1:n)*coefsz_r; 

f_L=leftSurface.faces;
v_L=new_leftSurface.vertices;

f_R=rightSurface.faces;
v_R=new_rightSurface.vertices;

figure;patch('Faces',f_L,'Vertices',v_L,'facecolor',[0.7,0.7,0.7],'edgecolor','none')
hold on
patch('Faces',f_R,'Vertices',v_R,'facecolor',[0.7,0.7,0.7],'edgecolor','none')
axis off; axis image; view([0 90]);set(gcf,'Color','w');
camlight
material dull
