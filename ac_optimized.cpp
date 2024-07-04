// compute the length of the edge ac under the assumption that the support sphere has radius r
// ac^2 is the root of a complicated polynomial DX^2+EX+F (where D<0)
// here, the formulas for D, E and F are optimized depending on the radii of the spheres
I ac(block B,bool verbose)
{
    I ab2=square(B[0]);
    I ad2=square(B[2]);
    I bc2=square(B[3]);
    I bd2=square(B[4]);
    I cd2=square(B[5]);
    //#if defined(uuuu)
    I D=square(bd2-I(4))-I(16);
    I E=-I(8)*square(bd2) + I(8)*(ad2 - I(4))*bc2 - I(2)*((ad2 - I(4))*bc2 - I(4)*ad2 - I(16))*bd2 - I(8)*(ad2 - I(4))*cd2;
    I F=(square(ad2) - I(8)*ad2)*square(bc2) - I(16)*square(cd2) - I(8)*(square(ad2) - I(4)*ad2)*bc2 + I(8)*(ad2*bc2 - I(4)*ad2)*bd2 - I(8)*((bc2 - I(4))*bd2 - I(4)*ad2 - I(4)*bc2 + I(16))*cd2;
    //#endif

    // ac^2 is a root of D*X^2+E*X+F
    if(verbose)
    {
    printf("D=%f,%f\n",lower(D),upper(D));
    printf("E=%f,%f\n",lower(E),upper(E));
    printf("F=%f,%f\n",lower(F),upper(F));
    }

    // possible length for ac
    I H=hull(ra+rc,ra+rc+I(2)*r);
    // no chance for Dx^2+Ex+F to have a root in H -> return empty
    if (!zero_in(D*pow(H,4)+E*pow(H,2)+F)) {//printf("no way\n");
	return I(2,1);}
    // negative discriminant: no such FM-tetrahedron
    I disc2=square(E)-I(4)*D*F;
    if (upper(disc2)<0) {//printf("discri<0\n");
	return I(2,1);}
    disc2=intersect(disc2,hull(I(0),I(INFINITY)));

//    if (zero_in(D)) printf("bla\n");

    // otherwise two potential solutions:
    I ac1=intersect(sqrt((-E+sqrt(disc2))/I(2)/D),H);
    I ac2=intersect(sqrt((-E-sqrt(disc2))/I(2)/D),H);
    if(empty(ac1)&&!empty(ac2)) {
    I ac3=sqrt((-E-sqrt(disc2))/I(2)/D);
//    printf("D=%f,%f\n",lower(D),upper(D));
//    printf("ac2=%f,%f\n",lower(ac3),upper(ac3));
    return ac2;}
    else if(empty(ac2)&&!empty(ac1)) {printf("ac1\n"); return ac1;}
    else if(empty(ac1)&&empty(ac2)) {printf("neither ac1 nor ac2\n"); return I(2,1);} // no such FM-tetrahedron
    else
    {
//        printf("ac1=%f,%f, ac2=%f,%f\n",lower(ac1),upper(ac1),lower(ac2),upper(ac2));
        return hull(ac1,ac2); // two solutions?
    }
    
    // if D contains 0 or is very small: use the continuity of the roots in the coefficients!
/*    if (zero_in(D))
    {
        printf("zero in D");
        I x=I(4)*D*F/square(E);
        x=hull(x,I(0)); // the interval x must contain 0
        if (upper(x)<0.78) =intersect(ac,-F/E*(I(1)+x));
    }
    */
}


