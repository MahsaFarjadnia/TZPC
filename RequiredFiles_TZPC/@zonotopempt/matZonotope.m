function [res] = matZonotope(obj,numofcol)
%MATZONOTOPE Summary of this function goes here
%   Detailed explanation goes here
newcen = vec2mat(obj.Z(:,1),numofcol);
for i =1:size(obj.generators,2)
    newgen{i} = vec2mat(obj.Z(:,i+1),numofcol);
end
res = matZonotope(newcen,newgen);
end