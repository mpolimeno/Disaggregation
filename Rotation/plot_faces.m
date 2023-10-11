function plot_faces(finalposint, finalndir,Nf,fignum,col)
% finalposint is the position of the center of each faces
% finalndir is the index of the coordinate normal to the faces
% finalori is useless
% Nf is the number of faces
% fignum is the number of the figure where you want to display the results.

figure(fignum)

for i=1:Nf

    if (finalndir(i) == 1)
        xs = [finalposint(i,1),   finalposint(i,1),   finalposint(i,1),   finalposint(i,1)];
        ys = [finalposint(i,2)+1, finalposint(i,2)+1, finalposint(i,2)-1, finalposint(i,2)-1];
        zs = [finalposint(i,3)+1, finalposint(i,3)-1, finalposint(i,3)-1, finalposint(i,3)+1];
        hold on
    end
    if (finalndir(i) == 2)
        xs = [finalposint(i,1)+1, finalposint(i,1)+1, finalposint(i,1)-1, finalposint(i,1)-1];
        ys = [finalposint(i,2),   finalposint(i,2),   finalposint(i,2),   finalposint(i,2)];
        zs = [finalposint(i,3)+1, finalposint(i,3)-1, finalposint(i,3)-1, finalposint(i,3)+1];
        hold on
    end
    if (finalndir(i) == 3)
        xs = [finalposint(i,1)+1, finalposint(i,1)-1, finalposint(i,1)-1, finalposint(i,1)+1];
        ys = [finalposint(i,2)+1, finalposint(i,2)+1, finalposint(i,2)-1, finalposint(i,2)-1];
        zs = [finalposint(i,3),   finalposint(i,3),   finalposint(i,3),   finalposint(i,3)];
        hold on
    end
    
    patch(xs,ys,zs,col)
end
xlabel('x',"Interpreter","LaTex")
ylabel('y',"Interpreter","LaTex")
zlabel('z',"Interpreter","LaTex")
set(gca,"FontSize",25)
axis equal
az = 121;
el = 31;
view(az,el)

end