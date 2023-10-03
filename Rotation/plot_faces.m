function y = plot_faces(finalposint, finalndir, finalori,Nf,fignum,col);
% finalposint is the position of the center of each faces
% finalndir is the index of the coordinate normal to the faces
% finalori is useless
% Nf is the number of faces
% fignum is the number of the figure where you want to display the results.

figure(fignum)
% clf

for i=1:Nf

	if (finalndir(i) == 1)
		xs = [finalposint(i,1) finalposint(i,1) finalposint(i,1) finalposint(i,1) ];
		ys =  [finalposint(i,2)+1 finalposint(i,2)+1 finalposint(i,2)-1 finalposint(i,2)-1 ];
		zs =  [finalposint(i,3)+1 finalposint(i,3)-1 finalposint(i,3)-1 finalposint(i,3)+1 ];
	%	plot3([finalposint(i,1) finalposint(i,1)+finalori(i)],[finalposint(i,2) finalposint(i,2)],[finalposint(i,3) finalposint(i,3)])
		hold on
	end;
	if (finalndir(i) == 2)
		xs = [finalposint(i,1)+1 finalposint(i,1)+1 finalposint(i,1)-1 finalposint(i,1)-1 ];
		ys =  [finalposint(i,2) finalposint(i,2) finalposint(i,2) finalposint(i,2) ];
		zs =  [finalposint(i,3)+1 finalposint(i,3)-1 finalposint(i,3)-1 finalposint(i,3)+1 ];
	%	plot3([finalposint(i,1) finalposint(i,1)],[finalposint(i,2) finalposint(i,2)+finalori(i)],[finalposint(i,3) finalposint(i,3)])
		hold on
	end;
	if (finalndir(i) == 3)
		xs = [finalposint(i,1)+1 finalposint(i,1)-1 finalposint(i,1)-1 finalposint(i,1)+1 ];
		ys =  [finalposint(i,2)+1 finalposint(i,2)+1 finalposint(i,2)-1 finalposint(i,2)-1 ];
		zs =  [finalposint(i,3) finalposint(i,3) finalposint(i,3) finalposint(i,3) ];
	%	plot3([finalposint(i,1) finalposint(i,1)],[finalposint(i,2) finalposint(i,2)],[finalposint(i,3) finalposint(i,3)+finalori(i)])
		hold on
	end;

%	patch(xs,ys,zs,'c')
%	patch(xs,ys,zs,[0.65 0.65 0.65])
		patch(xs,ys,zs,col)
end;

axis equal
%view(45,45)
y=0;

% tests
%np = 187;
%xs = [xc(np,1), xc(np,1),xc(np,1)+0,xc(np,1)+0];
%ys = [xc(np,2)+1, xc(np,2)+1,xc(np,2)-1,xc(np,2)-1];
%zs = [xc(np,3)+1, xc(np,3)-1,xc(np,3)-1,xc(np,3)+1];
%patch(xs,ys,zs,'r')