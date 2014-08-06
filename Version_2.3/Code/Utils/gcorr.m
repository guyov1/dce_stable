function Out=gcorr(A,B)
C=corrcoef(A,B);
Out=C(2);