function R = rot(fi, ax)
% This functions outputs a rotation matrix R. 
% The angle is given by "fi".
% The axis of rotation is given by "ax", 1 for X; 2 for Y and 3 for Z;
switch ax
    case 1
        R=[1 0 0;...
           0 cos(fi) sin(fi);...
           0 -sin(fi) cos(fi)];
    case 2
        R=[cos(fi) 0 -sin(fi);...
           0       1  0;...
           sin(fi) 0 cos(fi)];
    case 3
        R=[cos(fi) sin(fi) 0;...
           -sin(fi) cos(fi) 0;...
           0 0 1];
end

if (ax~=1)&& (ax~=2)&&(ax~=3)==1
    error('wrong axis!')
end

