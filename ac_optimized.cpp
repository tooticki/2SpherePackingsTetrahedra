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
    I D = (square(bd)*(I(-10)))*square(r)+pow(bd,4)+(pow(r,3)*(I(-12)))*rd+r*I(4)*pow(rd,3)+(square(r)*(I(-2)))*square(rd)+pow(r,4)*I(9)+(square(bd)*(I(-2)))*square(rd)+((square(bd)*(I(-4)))*r)*rd+pow(rd,4);
    I E =  ((square(r)*(I(-12)))*square(bc))*square(rd)+((pow(r,3)*(I(-8)))*square(bc))*rd+((r*(I(-8)))*square(bc))*pow(rd,3)+((square(ad)*(I(-2)))*square(bd))*square(bc)+square(cd)*I(40)*pow(r,4)+(pow(bd,4)*(I(-16)))*square(r)+square(bd)*I(40)*pow(r,4)+pow(r,4)*I(16)*square(rd)+pow(r,5)*I(96)*rd+(pow(r,3)*(I(-32)))*pow(rd,3)+(square(r)*(I(-8)))*pow(rd,4)+pow(r,6)*(I(-72))+(square(bc)*(I(-2)))*pow(rd,4)+square(bd)*I(2)*square(bc)*square(rd)+square(ad)*I(2)*square(bc)*square(rd)+square(cd)*I(16)*pow(r,3)*rd+square(cd)*I(8)*square(r)*square(rd)+(pow(r,4)*(I(-34)))*square(bc)+square(bd)*I(8)*square(r)*square(rd)+square(bd)*I(16)*pow(r,3)*rd+square(ad)*I(16)*square(bd)*square(r)+square(ad)*I(10)*square(r)*square(bc)+square(bd)*I(10)*square(r)*square(bc)+((square(ad)*(I(-16)))*square(cd))*square(r)+square(bd)*I(4)*r*square(bc)*rd+square(bd)*I(8)*square(cd)*square(r)+square(ad)*I(4)*r*square(bc)*rd;
    I F = pow(r,5)*I(96)*square(bc)*rd+((square(r)*(I(-8)))*square(bc))*pow(rd,4)+r*I(4)*pow(bc,4)*pow(rd,3)+((square(r)*(I(-2)))*pow(bc,4))*square(rd)+((pow(r,3)*(I(-12)))*pow(bc,4))*rd+pow(r,4)*I(16)*square(bc)*square(rd)+((pow(r,3)*(I(-32)))*square(bc))*pow(rd,3)+(((square(bd)*(I(-16)))*square(cd))*square(r))*square(bc)+(pow(cd,4)*(I(-48)))*pow(r,4)+(square(cd)*(I(-160)))*pow(r,6)+square(ad)*I(8)*square(cd)*square(r)*square(bc)+pow(r,5)*I(64)*pow(rd,3)+(pow(r,6)*(I(-32)))*square(rd)+(pow(r,7)*(I(-192)))*rd+pow(r,4)*I(16)*pow(rd,4)+pow(r,8)*I(144)+square(ad)*I(16)*square(bd)*square(r)*square(bc)+pow(rd,4)*pow(bc,4)+((square(ad)*(I(-2)))*pow(bc,4))*square(rd)+((square(cd)*(I(-32)))*pow(r,4))*square(rd)+((square(cd)*(I(-64)))*pow(r,5))*rd+pow(r,4)*I(9)*pow(bc,4)+(pow(r,6)*(I(-72)))*square(bc)+((square(ad)*(I(-64)))*square(bd))*pow(r,4)+square(cd)*I(8)*square(r)*square(bc)*square(rd)+square(cd)*I(16)*pow(r,3)*square(bc)*rd+((square(ad)*(I(-10)))*square(r))*pow(bc,4)+((pow(ad,4)*(I(-16)))*square(r))*square(bc)+square(ad)*I(40)*pow(r,4)*square(bc)+square(ad)*I(64)*square(cd)*pow(r,4)+pow(bc,4)*pow(ad,4)+square(cd)*I(40)*pow(r,4)*square(bc)+square(bd)*I(64)*square(cd)*pow(r,4)+square(ad)*I(16)*pow(r,3)*square(bc)*rd+square(ad)*I(8)*square(r)*square(bc)*square(rd)+(((square(ad)*(I(-4)))*r)*pow(bc,4))*rd;
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
