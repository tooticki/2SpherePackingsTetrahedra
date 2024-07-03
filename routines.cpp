// basic functions

// return false if x<y+z (or one of the permutation) does not hold 
bool trig(I x, I y, I z)
{
    return (lower(x)<upper(y+z) && lower(y)<upper(x+z) && lower(z)<upper(x+y));
}

// return false if the facial condition is not satisfied
bool facial(I ab, I ac, I ad, I bc, I bd, I cd)
{
    return (trig(ab,ac,bc) && trig(ab,bd,ad) && trig(ac,cd,ad) && trig(bc,bd,cd));
}

// angle in A in triange ABC (law of cosines)
I angle(I bc, I ab, I ac)
{
    return acos(intersect((square(ab)+square(ac)-square(bc))/(I(2)*ac*ab),I(-1,1)));
}

// dihedral angle along ab (spherical law of cosines)
I dihedral(I ab, I ac, I ad, I bc, I bd, I cd)
{
    I ab2=square(ab);
    I ac2=square(ac);
    I ad2=square(ad);
    I bc2=square(bc);
    I bd2=square(bd);
    I cd2=square(cd);
    return acos(((ac2+ad2+bc2+bd2-I(2)*cd2-ab2)*ab2 + (bc2-ac2)*(ad2-bd2))/sqrt((square(ab+ac)-bc2)*(square(ab-ac)-bc2)*(square(ab+ad)-bd2)*(square(ab-ad)-bd2)));
}

// return false if the cayley-menger nonnegative determinant is not satisfied
bool cayley_menger(I ab, I ac, I ad, I bc, I bd, I cd)
{
    return (upper(square(ab*ac*bd)+square(ab*ac*cd)+square(ab*ad*bc)+square(ab*ad*cd)+square(ab*bc*cd)+square(ab*bd*cd)+square(ac*ad*bc)+square(ac*ad*bd)+square(ac*bc*bd)+square(ac*bd*cd)+square(ad*bc*bd)+square(ad*bc*cd)-square(ab*ab*cd)-square(ab*ac*bc)-square(ab*ad*bd)-square(ab*cd*cd)-square(ac*ac*bd)-square(ac*ad*cd)-square(ac*bd*bd)-square(ad*ad*bc)-square(ad*bc*bc)-square(bc*bd*cd))>0);
}

// return false if edge lengths are not compatible with a tetrahedron
bool tetrahedral(I ab, I ac, I ad, I bc, I bd, I cd)
{
    return (facial(ab,ac,ad,bc,bd,cd) && cayley_menger(ab,ac,ad,bc,bd,cd));
}


// volume of the tetrahedron
I vol(I ab, I ac, I ad, I bc, I bd, I cd)
{
    I ab2=square(ab);
    I ac2=square(ac);
    I ad2=square(ad);
    I bc2=square(bc);
    I bd2=square(bd);
    I cd2=square(cd);
    return sqrt((I(4)*ad2*ab2-square(ab2+ad2-bd2))*(I(4)*ac2*ab2-square(ab2+ac2-bc2))-square((ac2+ad2-cd2)*I(2)*ab2-(ab2+ac2-bc2)*(ab2+ad2-bd2)))/I(24)/ab; // general case: 21 variables
}

/***************************************************/

// solid angle in A : Lagrange formula
I solid1(I ab, I ac, I ad, I bc, I bd, I cd)
{
    I ab2=square(ab);
    I ac2=square(ac);
    I ad2=square(ad);
    I bc2=square(bc);
    I bd2=square(bd);
    I cd2=square(cd);
    I V=vol(ab,ac,ad,bc,bd,cd);
    I f=(I(2)*ab*ac*ad+(ac2+ad2-cd2)*ab+(ab2+ad2-bd2)*ac+(ab2+ac2-bc2)*ad)/I(12);
    I z=I(2)*atan(V/f);
    if (lower(z)>=0) return z;
    if (upper(z)<=0) return z+pi_twice<I>();
    // otherwise unclear...solid angle in two disjoint intervals [0,a] and [b,pi] ?
//    return intersect(z,I(0,4)); // discard [b,pi] (why?)
    return hull(I(0),pi<I>());
}

// solid angle in A : Girard
I solid2(I ab, I ac, I ad, I bc, I bd, I cd)
{
    return dihedral(ab,ac,ad,bc,bd,cd)+dihedral(ac,ad,ab,cd,bc,bd)+dihedral(ad,ab,ac,bd,cd,bc)-pi<I>();
}

// solid angle : LHuiler formula
I solid3(I ab, I ac, I ad, I bc, I bd, I cd)
{
    I t=(angle(bc,ab,ac)+angle(cd,ac,ad)+angle(bd,ab,ad))/I(2);
    I z=tan(t/I(2))*tan((t-angle(bc,ab,ac))/I(2))*tan((t-angle(cd,ac,ad))/I(2))*tan((t-angle(bd,ab,ad))/I(2));
    z=sqrt(intersect(z,I(0,INFINITY))); // better way to do ?
    z=I(4)*atan(z);
    if (lower(z)>=0) return z;
    if (upper(z)<=0) return z+I(4)*pi<I>();
    return hull(I(0),pi<I>());
}

// combine the three formulas
I solid(I ab, I ac, I ad, I bc, I bd, I cd)
{
    I z=solid1(ab,ac,ad,bc,bd,cd);
    z=intersect(z,solid2(ab,ac,ad,bc,bd,cd));
    z=intersect(z,solid3(ab,ac,ad,bc,bd,cd));
    return z;
}

/***************************************************/


// volume of the part covered by balls
I cov(I ab, I ac, I ad, I bc, I bd, I cd)
{
    return
    solid(ab,ac,ad,bc,bd,cd)/I(3)*pow(ra,3)+
    solid(bc,bd,ab,cd,ac,ad)/I(3)*pow(rb,3)+
    solid(cd,ac,bc,ad,bd,ab)/I(3)*pow(rc,3)+
    solid(ad,bd,cd,ab,ac,bc)/I(3)*pow(rd,3);
}

/***************************************************/


// density
I density1(I ab, I ac, I ad, I bc, I bd, I cd)
{
    return cov(ab,ac,ad,bc,bd,cd)/vol(ab,ac,ad,bc,bd,cd);
}

// alternative formula relying on atan(x)<=x for x>=0, useful for small volumes
I aux(I ab, I ac, I ad, I bc, I bd, I cd)
{
    I ab2=square(ab);
    I ac2=square(ac);
    I ad2=square(ad);
    I bc2=square(bc);
    I bd2=square(bd);
    I cd2=square(cd);
    return (I(2)*ab*ac*ad+(ac2+ad2-cd2)*ab+(ab2+ad2-bd2)*ac+(ab2+ac2-bc2)*ad)/I(12);
}
I density2(I ab, I ac, I ad, I bc, I bd, I cd)
{
    I V=vol(ab,ac,ad,bc,bd,cd);
    I d=0;
    I f;
    f=aux(ab,ac,ad,bc,bd,cd);
    if (lower(V/f)>=0) d+=I(2)*hull(I(0),I(1)/f)*pow(ra,3)/I(3);
    else d+=I(2)*(atan(V/f)+pi<I>())/V*pow(ra,3)/I(3);
    f=aux(bc,bd,ab,cd,ac,ad);
    if (lower(V/f)>=0) d+=I(2)*hull(I(0),I(1)/f)*pow(rb,3)/I(3);
    else d+=I(2)*(atan(V/f)+pi<I>())/V*pow(rb,3)/I(3);
    f=aux(cd,ac,bc,ad,bd,ab);
    if (lower(V/f)>=0) d+=I(2)*hull(I(0),I(1)/f)*pow(rc,3)/I(3);
    else d+=I(2)*(atan(V/f)+pi<I>())/V*pow(rc,3)/I(3);
    f=aux(ad,bd,cd,ab,ac,bc);
    if (lower(V/f)>=0) d+=I(2)*hull(I(0),I(1)/f)*pow(rd,3)/I(3);
    else d+=I(2)*(atan(V/f)+pi<I>())/V*pow(rd,3)/I(3);
    return d;
}

// combine both formulas
I density(I ab, I ac, I ad, I bc, I bd, I cd)
{
//    return density1(ab,ac,ad,bc,bd,cd);
    return intersect(density1(ab,ac,ad,bc,bd,cd),density2(ab,ac,ad,bc,bd,cd));
}


