function m=normalized_Mask(var,i,j)
m=-var+(var-(-var)).*rand(i,j);
m=m./norm(m);
end
