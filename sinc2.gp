\\ Example Calculation by Bradley Klee 
\\ Definite integral of cardinal sine squared
\\ 	Int_{0..10000*Pi}(sinc(x))^2*dx
\\ 	where sinc(x) = sin(x)/x  

\\ input: integer for series truncation
\\ output: expansion coefficients of (sinc(x))^2 around x=0
{sinc20(mPrec)=my(cv);cv=[1.,0.];for(j=2,mPrec,
	cv=concat(cv,-[cv[j-1]*4/(2+3*j+j^2)]));
	cv};

\\ input: integers for series truncation and expansion center
\\ output: expansion coefficients of (sinc(x))^2 around x=n*pi
{sinc2n(mPrec,nDom)=my(cv);
	cv=[0.,0.,1/(nDom*Pi)^2,-2/(nDom*Pi)^3];
	for(j=4,mPrec,cv=concat(cv,cv[-4..-1]
		*[4/((j-j^2)*nDom^2*Pi^2), 
		8/((1-j)*j*nDom*Pi), 
		(j-j^2-4*nDom^2*Pi^2)/((-j+j^2)*nDom^2*Pi^2), 
		-2/(nDom*Pi)]~));
	cv};

\\ input: integer for series truncation
\\ output: integrated powers of x, evaluated at -Pi
{xPiInt(mPrec)=vector(mPrec+1,j,(-Pi)^j/j)~};

\\ input: integers for total domain and series truncation	
\\ output: Int_{0..n*Pi}(sinc(x))^2*dx, with error estimate
\\ comment: this function is not completely optimized!
\\		Running sums would be better. 
{sinc2Int(nDom,mPrec)=my(IntVals,xVals);
	xVals=xPiInt(mPrec);
	IntVals=[
		sinc20(mPrec)*abs(xVals),
		-sinc2n(mPrec,1)*xVals];
	for(j=1,nDom-1,IntVals=matconcat([IntVals;[
		sinc2n(mPrec,j)*abs(xVals),
		-sinc2n(mPrec,j+1)*xVals]]));
	vector(nDom,j,1)*IntVals*[1/2,1;1/2,-1]};

\\ input: integer for total domain	
\\ output: Int_{0..n*Pi}(sinc(x))^2*dx
{sinc2IntIdiom(nDom)=
	intnum(x=0,nDom*Pi,(sin(x)/x)^2)};

\\ input: integer for total domain	
\\ output: Int_{0..n*Pi}(sinc(x))^2*dx
{sinc2IntPiecewiseIdiom(nDom)=my(int);int=0;
	for(n=0,nDom-1,int+=intnum(x=n*Pi,(n+1)*Pi,
		(sin(x)/x)^2));
	int};

\p 50
Pi/2
sinc2IntIdiom(10000) /* not accurate */
sinc2IntPiecewiseIdiom(10000)
out=sinc2Int(10000,75)
log(abs(Pi/2-out[1]))

