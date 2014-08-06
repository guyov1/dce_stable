function Out=gpowspace(Mn,Mx,n,Pow)
Out=Mn+(Mx-Mn)*(linspace(0,1,n).^Pow);