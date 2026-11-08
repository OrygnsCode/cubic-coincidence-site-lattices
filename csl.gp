\\ Number of coincidence site lattices of index 2n-1 in Z^3, up to the rotations of the cube.
\\ Closed formula from P. Zeiner, Z. Kristallogr. 220 (2005), Theorem 2.2.
\\ Usage in gp:  \r csl.gp   then   vector(30, n, a(n))   or   writeterms(10000)

a(n) = {my(S = 2*n-1, f, p, k, psi, n1, n2 = 0, n3 = 0, n4 = 0, big);
  if(S == 1, return(1));
  f = factor(S); p = f[,1]; k = #p; psi = S*prod(i=1, k, 1+1/p[i]);
  n1 = (S == 3);
  if(prod(i=1, k, p[i]%4 == 1), n2 = 2^(k-1));
  big = select(q -> q != 3, p~);
  if(S > 3 && #big && prod(i=1, #big, big[i]%6 == 1) && valuation(S, 3) <= 1, n3 = 2^(#big-1));
  if(S > 3 && prod(i=1, k, p[i]%8 == 1 || p[i]%8 == 3), n4 = 2^(k-1));
  n1 + n2 + n3 + n4 + (psi - 4*n1 - 6*n2 - 8*n3 - 12*n4)/24;}

writeterms(N) = for(n=1, N, write("terms.txt", n, " ", a(n)));
