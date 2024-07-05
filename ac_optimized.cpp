// compute the length of the edge ac under the assumption that the support sphere has radius r
// ac^2 is the root of a complicated polynomial DX^2+EX+F (where D<0)
// here, the formulas for D, E and F are optimized depending on the radii of the spheres
I ac(block B, bool verbose=false)
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

    // no chance for Dx^2+Ex+F to have a root in H^2 -> return empty
    if (!zero_in(D*pow(H,4)+E*pow(H,2)+F)) {return I(1,0);}

    // usual formula through discriminant
    I disc2=square(E)-I(4)*D*F;
    if (upper(disc2)<=0) return I(1,0); // no FM-tetrahedron

    // otherwise two potential solutions for ac^2, the good one seems to always be this one
    I acs2=(-E-sqrt(disc2))/I(2)/D;


    // if D contains 0 or is very small: use the continuity of the roots in the coefficients!
    if (zero_in(D)||zero_in(F))
	{
	    I x=I(4)*D*F/square(E);
	    x=hull(x,I(0)); // the interval x must contain 0
	    if (upper(x)<0.78)
		{
		    I acs3=intersect(acs2,-F/E*(I(1)+x));
		    if (width(acs3)<width(acs2)) {acs2=acs3;}
		}
	}


    I acs=sqrt(acs2);
    acs=intersect(acs,H);
    return acs;
}
