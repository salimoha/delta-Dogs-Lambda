function [ind]=Search_triangulation_neighber_find(ind,tri,xi,yil);
% The code to find the negibers of a point with index ind

%tri=delaunayn(xi.');
global n

[row_indices,column_indices]=find(tri==ind);
[row_indices1,column_indices1]=find(tri<n+2);
row_indices=setdiff(row_indices,row_indices1);
xil=xi(:,unique(tri(row_indices,:)));
yil=yi(unique(tri(row_indices,:)));



end