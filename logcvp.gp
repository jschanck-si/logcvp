/*
 * Proof of Concept implementation of the CVP subroutine
 * for finding a short generator of a principal ideal in
 * prime and power-of-two cyclotomic fields.
 *
 * usage from gp:
 *        read("logcvp.gp")
 *        test([dimension]) \\ Default 512
 * The first test will perform some precomputation that is
 * reused in subsequent tests at the same dimension.
 *
 * author: John M. Schanck <jschanck@uwaterloo.ca>
 * date: 2015-03-09
 * license: CC0 <https://creativecommons.org/publicdomain/zero/1.0/>
 * repo: https://github.com/jschanck-si/logcvp
 */


\p 150

proj(u, v) = {
  (u~ * v)/(u~ * u) * u;
}

\\ Gram-Schmidt orthogonalization
gs(M) = {
  my(M2);
  M2 = matconcat([M[,1], M[,2] - proj(M[,1], M[,2])]);
  for(j=3, #M[1,],
    M2 = concat(M2, M[,j] - vecsum(vector(j-1, i, proj(M2[,i], M[,j]))))
  );
  M2;
}

\\ Plot log base q lengths of M's GSO vectors
basisProfile(M,q=2) = {
  my(GM, V);
  GM = gs(1.*M);
  V = log(vector(#M, i, sqrt(norml2(GM[,i]))))/log(q);
  plothraw(vector(#V,i,i), V, 1);
  V;
}

\\ Babai's nearest planes algorithm
babai(y, T, G) = {
  my(N, x, u, k);
  N = #y~;
  x = y;
  u = vector(N)~;
  for(i=0, N-1,
      j=N-i;
      u[j] = round((x~ * G[,j])/(norml2(G[,j])));
      x = x - u[j]*T[,j]);
  [x, u];
}

\\ Sample a discrete gaussian on ZZ centered at c with std dev s
sampleZ(n, s, c) = {
  my(t,d,low,high,q,r);
  t = log(n);
  d = -Pi/s^2;
  low = floor(c - t*s);
  high = ceil(c + t*s);
  while(1,
    q = low + random(high - low + 1);
    r = exp(d*(q-c)^2);
    if(random(1.) <= r, break));
  q;
}

\\ Sample polynomial of degree n with coefficients from sampleZ(n,s,0)
sampleGaussR(n, s) = {
  Pol(vector(n,i,sampleZ(n,s,0)));
}

\\ Log embedding for totally complex fields
embed(v, roots) = {
  vector(#roots, j, 2*log(abs(subst(lift(v),x,roots[j]))))~;
}

cvpPrep(N) = {
  my(r2, U, roots, M, M2, T, GM);
  if(N%2 == 0 && isprimepower(N), \\ N is a power of 2
    r2 = eulerphi(N)/2;
    U = vector(r2-1, i, Mod(x^((1-(2*i+1))/2)*(1 - x^(2*i+1))/(1-x), polcyclo(N)));
    roots = vector(r2, j, exp(Pi*I*(2*j-1)/(2*r2)));
    M = matconcat(vector(#U, i, embed(U[i],roots))),
  if(isprime(N), \\ N is prime
    r2 = eulerphi(N)/2;
    U = vector(r2-1, i, Mod(x^lift(Mod((1-(i+1))/2,N))*(1 - x^(i+1))/(1-x), polcyclo(N)));
    roots = vector(r2, j, exp(2*Pi*I*j/N));
    M = matconcat(vector(#U, i, embed(U[i],roots))),
  error("Cannot handle N=", N)));
  \\ Explicit unit group computation, could be useful for small composite N
  \\M = matconcat(vector(#K.fu, i, embed(K.fu[i],roots)));

  \\TODO: try to determine whether cyclotomic units are index 1
  \\Reg = abs(matdet(matrix(r2-1,r2-1,i,j,M[i,j])));

  \\ Preprocess basis with LLL
  M = matrix(r2-1, r2-1, i,j,M[i,j]);
  T = qflll(M);
  M = M * T;
  GM = gs(M);

  \\ No preprocessing
  \\M = matrix(r2-1, r2-1, i,j,M[i,j]);
  \\GM = gs(M);
  \\T = matid(#M);

  [N, r2, M, GM, roots, U, T];
}

\\ Sample unit with using Gaussian distribution on exponent vector
randUnit(prep,s=3) = {
  my(R,RU,N,r2,M,GM,roots,U,T);
  [N, r2, M, GM, roots, U, T] = prep;
  R = vector(#U,i,sampleZ(#U,s,0));
  RU = lift(Mod(prod(i=1, #R, U[i]^R[i]), polcyclo(N)));
}

\\ Attempt CVP in Log unit lattice for target Log(f) - vecavg(Log(f))
\\ flags: 1 - print recovered exponent vector
smallGenerator(f, prep, flags=0)={
  my(c1,o,u,N, r2, M,GM,roots,U,T);
  [N, r2, M, GM, roots, U, T] = prep;

  c1 = embed(f, roots);
  c1 = vector(r2-1, i, c1[i]);
  c1 -= vecsum(c1)*vector(#c1,i,1/(r2));

  [o, u] = babai(c1~, M, GM);
  u = T*u;
  if(bitand(flags,1),printf("%s\n",u));

  f*prod(i=1,#u,U[i]^(-u[i]));
}

\\ flags: 1 - print exponent vectors for smallGenerator(r), smallGenerator(f), and smallGenerator(g)
test(N=512,flags=0) = {
  \\ Call cvpPrep, if prep does not exist or prep[1] != N
  iferr([Np,r2,M,GM,roots,U,T] = prep, E, prep=cvpPrep(N); [Np,r2,M,GM,roots,U,T] = prep);
  if(N != Np, prep=cvpPrep(N); [Np,r2,M,GM,roots,U,T] = prep);
  m = eulerphi(N);

  s = sqrt(128*m);
  f = sampleGaussR(m, s);
  \\ Uncomment to ensure f has prime norm
  \\until(isprime(polresultant(f, polcyclo(N))), f = sampleGaussR(m, s));
  f = Mod(f, polcyclo(N));

  r = randUnit(prep,4);
  fr = f*r;
  \\ Uncomment here and in smallGenerator to see exponent vectors of r, f and g
  if(bitand(flags,1),
    printf("r ~ ");smallGenerator(r,prep,flags);
    printf("f ~ ");smallGenerator(f,prep,flags);
    printf("g ~ "));
  g = smallGenerator(fr,prep,flags);
  printf("\nnorm initial generator = %d\nnorm random generator = %d\nnorm recovered generator: %d\n",\
    round(sqrt(norml2(lift(f)))),\
    round(sqrt(norml2(lift(fr)))),\
    round(sqrt(norml2(lift(g)))));
  printf("Difference (|found| - |orig|) %d", round(sqrt(norml2(lift(g)))) - round(sqrt(norml2(lift(f)))));
}

