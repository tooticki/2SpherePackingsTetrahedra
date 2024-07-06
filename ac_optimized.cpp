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
#if defined(uuuu) && defined(support_sphere_r)
    I D=square(bd2-I(4))-I(16);
    I E=-I(8)*square(bd2) + I(8)*(ad2 - I(4))*bc2 - I(2)*((ad2 - I(4))*bc2 - I(4)*ad2 - I(16))*bd2 - I(8)*(ad2 - I(4))*cd2;
    I F=(square(ad2) - I(8)*ad2)*square(bc2) - I(16)*square(cd2) - I(8)*(square(ad2) - I(4)*ad2)*bc2 + I(8)*(ad2*bc2 - I(4)*ad2)*bd2 - I(8)*((bc2 - I(4))*bd2 - I(4)*ad2 - I(4)*bc2 + I(16))*cd2;
#elif defined(rrrr) && defined(support_sphere_r) // TODO simplify formulas?
    I D = (bd2*(I(-10)))*square(r)+pow(bd2,2)+(pow(r,3)*(I(-12)))*rd+r*I(4)*pow(rd,3)+(square(r)*(I(-2)))*square(rd)+pow(r,4)*I(9)+(bd2*(I(-2)))*square(rd)+((bd2*(I(-4)))*r)*rd+pow(rd,4);
    I E =  ((square(r)*(I(-12)))*bc2)*square(rd)+((pow(r,3)*(I(-8)))*bc2)*rd+((r*(I(-8)))*bc2)*pow(rd,3)+((ad2*(I(-2)))*bd2)*bc2+cd2*I(40)*pow(r,4)+(pow(bd2,2)*(I(-16)))*square(r)+bd2*I(40)*pow(r,4)+pow(r,4)*I(16)*square(rd)+pow(r,5)*I(96)*rd+(pow(r,3)*(I(-32)))*pow(rd,3)+(square(r)*(I(-8)))*pow(rd,4)+pow(r,6)*(I(-72))+(bc2*(I(-2)))*pow(rd,4)+bd2*I(2)*bc2*square(rd)+ad2*I(2)*bc2*square(rd)+cd2*I(16)*pow(r,3)*rd+cd2*I(8)*square(r)*square(rd)+(pow(r,4)*(I(-34)))*bc2+bd2*I(8)*square(r)*square(rd)+bd2*I(16)*pow(r,3)*rd+ad2*I(16)*bd2*square(r)+ad2*I(10)*square(r)*bc2+bd2*I(10)*square(r)*bc2+((ad2*(I(-16)))*cd2)*square(r)+bd2*I(4)*r*bc2*rd+bd2*I(8)*cd2*square(r)+ad2*I(4)*r*bc2*rd;
    I F = pow(r,5)*I(96)*bc2*rd+((square(r)*(I(-8)))*bc2)*pow(rd,4)+r*I(4)*pow(bc2,2)*pow(rd,3)+((square(r)*(I(-2)))*pow(bc2,2))*square(rd)+((pow(r,3)*(I(-12)))*pow(bc2,2))*rd+pow(r,4)*I(16)*bc2*square(rd)+((pow(r,3)*(I(-32)))*bc2)*pow(rd,3)+(((bd2*(I(-16)))*cd2)*square(r))*bc2+(pow(cd2,2)*(I(-48)))*pow(r,4)+(cd2*(I(-160)))*pow(r,6)+ad2*I(8)*cd2*square(r)*bc2+pow(r,5)*I(64)*pow(rd,3)+(pow(r,6)*(I(-32)))*square(rd)+(pow(r,7)*(I(-192)))*rd+pow(r,4)*I(16)*pow(rd,4)+pow(r,8)*I(144)+ad2*I(16)*bd2*square(r)*bc2+pow(rd,4)*pow(bc2,2)+((ad2*(I(-2)))*pow(bc2,2))*square(rd)+((cd2*(I(-32)))*pow(r,4))*square(rd)+((cd2*(I(-64)))*pow(r,5))*rd+pow(r,4)*I(9)*pow(bc2,2)+(pow(r,6)*(I(-72)))*bc2+((ad2*(I(-64)))*bd2)*pow(r,4)+cd2*I(8)*square(r)*bc2*square(rd)+cd2*I(16)*pow(r,3)*bc2*rd+((ad2*(I(-10)))*square(r))*pow(bc2,2)+((pow(ad2,2)*(I(-16)))*square(r))*bc2+ad2*I(40)*pow(r,4)*bc2+ad2*I(64)*cd2*pow(r,4)+pow(bc2,2)*pow(ad2,2)+cd2*I(40)*pow(r,4)*bc2+bd2*I(64)*cd2*pow(r,4)+ad2*I(16)*pow(r,3)*bc2*rd+ad2*I(8)*square(r)*bc2*square(rd)+(((ad2*(I(-4)))*r)*pow(bc2,2))*rd;
#else//TODO other cases
    printf("Error: ac formula are not implemented for this case");
    I D = I(0);
    I E = I(0);
    I F = I(0);
#endif



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
