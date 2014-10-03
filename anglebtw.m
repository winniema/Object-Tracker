function[theta] = anglebtw(x1,y1,x2,y2)
% Input: (x1,y1) is the direction vector of the first line segment
% (x2,y2) is the direction vector of the second line segment

% Ouput: Angle. Returns the angle between the two line segments. The
% calculations is done using dot product

dotpro = (x1)*(x2)+(y1)*(y2);
mag_line1 = sqrt((x1)^2 + (y1)^2);
mag_line2 = sqrt((x2)^2 + (y2)^2);

cos_theta = dotpro/(mag_line1*mag_line2);
theta = acos(cos_theta);

end