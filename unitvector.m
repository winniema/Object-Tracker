function[i, j] = unitvector(x, y)
% Input: Direction vector
% Output: The Unit vector

    mag = sqrt(x^2+y^2);
    i = x/mag;
    j = y/mag;

end