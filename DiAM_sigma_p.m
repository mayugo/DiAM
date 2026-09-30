sigma =	[ 900 200 0; 
		    200 300 0; 
		        0  	0 0];               % Estat de tensions
[v,p] = eig(sigma);                 % Calcular valors 'p' i vectors propis 'v'
[P,ind] = sort(diag(p),'descend');  % Ordena valors de major a menor
p = p(ind,ind)
v = v(:,ind)