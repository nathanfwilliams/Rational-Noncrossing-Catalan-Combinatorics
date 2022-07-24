################################
# The Coxeter Group: Section 3 #
################################
W:=CoxeterGroup("I",2,8); # W is a finite irreducible Coxeter group
r:=W.nbGeneratingReflections; # r is the rank of W
h:=Order(W,EltWord( W , [1..r] )); # h is the Coxeter number
p:=1; # p is the parameter (coprime to h) in the formula Cat_p(W;q) := \prod_{i=1}^r \frac{[p+(p e_i \mod h)]_q}{[d_i]_q}

################################
# The Hecke Algebra: Section 6 #
################################
q:=X(Cyclotomics);;q.name:="q";;q:=Mvp(q);
H:=Hecke(W,q); # H is the Hecke algebra of W
Schurs:=SchurElements(H);

dim:=FakeDegrees(W,1); # the dimensions of the irreps of W
Feg:=FakeDegrees(W,q); # the fake degrees of the irreps of W (Equation 6.2)
FegZ:=FakeDegrees(W,E(h)^(p)); # the fake degrees of the irreps of W, evaluated at a primitive h-th root of unity; the fake degrees evaluated at a primitive h-th root of unity vanish except h many irreps (Remark ???)
Deg:=List(SchurElements(H),i->Schurs[PositionId(W)]/Mvp(i)); # the generic degrees of the irreps of H (Equation 6.3)
DegZ:=List(SchurElements(H),i->Value(Schurs[PositionId(W)]/Mvp(i),["q",E(h)^p])); # the generic degrees of the irreps of H, evaluated at a primitive h-th root of unity; the generic degrees evaluated at a primitive h-th root of unity vanish except on these r+1 exterior powers of the reflection representation (Remark 6.9)

ExtPowers:=ChevieCharInfo(W).extRefl; # the indices (among all irreps of W) of the exterior powers \wedge^i V (for i=0,1,\ldots,r) of the reflection representation V of W
Print("Remark 6.9: Fake degrees evaluated at a primitive h-th root of unity vanish except on the r+1 exterior powers of the reflection representations: ",Positions(List([1..Length(DegZ)],i->(i in ExtPowers)),true)=Positions(List([1..Length(DegZ)],i->DegZ[i]<>0),true),"\n");

# Compute the contents c(\chi) of the irreducibles \chi (Equation 6.6)
classes:=[];
for t in W.reflections do
	class:=ConjugacyClass(W,t);
	classpos:=Position(ConjugacyClasses(W),class);
	classsize:=Size(class);
	if not([classsize,classpos] in classes) then
		Add(classes,[classsize,classpos]);
	fi;
od;
Contents:=List(CharTable(W).irreducibles,chi->1/chi[1]*Sum(classes,i->i[1]*chi[i[2]]));

####################################################################
# Theorem 1.9.(1) and Corollary 6.11 (Catalan number computations) #
####################################################################
# Compute the rational q-Catalan number from Equation 1.3: Cat_p(W;q) := \prod_{i=1}^r \frac{[p+(p e_i \mod h)]_q}{[d_i]_q}
Catp:=Product([1..r],i->Sum([0..(p+Mod((W.degrees[i]-1)*p,h))-1],j->q^j))/(Schurs[PositionId(W)]);

# Compare Cat_p(W;q) with (1-q)^(-r) q^{rp/2}\sum_\chi q^{-p/h c(\chi)}\Feg_\chi (\zeta^p) \Deg_\chi (q), where \zeta is a primitive h-th root of unity
Print("The normalized LHS of the sum (\star) with Feg specialized to \zeta^p equals Cat_p(W;q): ",q^(r*p/2)/((1-q)^r*Schurs[PositionId(W)])*Sum([1..Length(Feg)],i->q^(-p/h*Contents[i])*FegZ[i]*Deg[i])=Catp,"\n");

# Compare Cat_p(W;q) with q^{rp/2}\sum_\chi q^{-p/h c(\chi)}\Feg_\chi (q) \Deg_\chi (\zeta^p), where \zeta is a primitive h-th root of unity
Print("The normalized RHS of the sum (\star) with Deg specialized to \zeta^p equals Cat_p(W;q): ",q^(r*p/2)/((1-q)^r*Schurs[PositionId(W)])*Sum([1..Length(Feg)],i->q^(-p/h*Contents[i])*Feg[i]*DegZ[i])=Catp,"\n");

#################################################################
# Theorem 1.9.(2) and Section 6.6 (Parking number computations) #
#################################################################
# Compute the rational q-Parking number from Equation 1.3: [p]_q^r
Parkp:=Sum([0..p-1],j->q^j)^r;

# Compare [p]_q^r with  (1-q)^(-r) q^{rp/2} \sum_\chi  q^{-p/h c(\chi)} \dim(\chi) \Feg_\chi (\zeta^p), where \zeta is a primitive h-th root of unity
Print("The normalized parking sum with Feg specialized to \zeta^p equals [p]^r: ",q^(r*p/2)/(1-q)^r*Sum([1..Length(Feg)],i->q^(-p/h*Contents[i])*FegZ[i]*dim[i])=Parkp,"\n");

# Compare [p]_q^r with (1-q)^(-r) q^{rp/2} \sum_\chi  q^{p/h c(\chi)} \dim(\chi) \Deg_\chi (\zeta^p), where \zeta is a primitive h-th root of unity
Print("The normalized parking sum with Deg specialized to \zeta^p equals [p]^r: ",q^(r*p/2)/(1-q)^r*Sum([1..Length(Feg)],i->q^(-p/h*Contents[i])*DegZ[i]*dim[i])=Parkp,"\n");
