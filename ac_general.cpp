// compute the length of the edge ac under the assumption that the support sphere has radius r
// ac^2 is the root of a complicated polynomial DX^2+EX+F (where D<0)
I ac(block B,bool verbose)
{
    I ab=B[0];
    I ad=B[2];
    I bc=B[3];
    I bd=B[4];
    I cd=B[5];
    I D=(((rd+bd-rb)*(rb+bd-rd))*(r*(I(-2))+bd-rb-rd))*(r*I(2)+bd+rb+rd);
    I E=rd*I(4)*(r*pow(ra,2)*(pow(bd,2)+pow(bc,2)-pow(cd,2))+ra*I(2)*pow(r,2)*(pow(bd,2)+pow(bc,2)-pow(cd,2))+(r*(pow(bd,2)+pow(ab,2)-pow(ad,2)))*pow(rc,2)+((pow(bd,2)*I(2)+pow(ab,2)*I(2)+pow(ad,2)*(I(-2)))*pow(r,2))*rc+(pow(bd,2)*(pow(bc,2)+pow(ab,2))+pow(cd,2)*pow(ab,2)-(pow(ab,2)*I(2)-pow(ad,2))*pow(bc,2))*r+((pow(bc,2)*I(2)+pow(ab,2)*I(2)+pow(bd,2)*(I(-4))+pow(ad,2)*I(2)+pow(cd,2)*I(2))*pow(r,2))*rb+(r*(pow(bc,2)+pow(ab,2)+pow(bd,2)*(I(-2))+pow(ad,2)+pow(cd,2)))*pow(rb,2))+(pow(cd,2)*(I(-2))+pow(ad,2)*(I(-2)))*pow(rb,4)+((pow(cd,2)*(I(-8))+pow(ad,2)*(I(-8)))*r)*pow(rb,3)+(((ra*(I(-16)))*pow(r,2))*pow(bd,2)+((pow(ra,2)*(I(-8)))*r)*pow(bd,2)+((r*(I(-4)))*pow(rb,2))*(pow(ab,2)-pow(bd,2)-pow(ad,2))+((pow(r,2)*(I(-8)))*rb)*(pow(ab,2)-pow(bd,2)-pow(ad,2))+(r*(I(-4)))*(pow(bd,4)-pow(bd,2)*(pow(ad,2)+pow(ab,2))))*rc+(((ra*(I(-8)))*r)*pow(bd,2)+(pow(rb,2)*(I(-2)))*(pow(ab,2)-pow(bd,2)-pow(ad,2))+(pow(ra,2)*(I(-4)))*pow(bd,2)+(pow(ad,2)*I(2)+pow(ab,2)*I(2))*pow(bd,2)+pow(bd,4)*(I(-2))+((r*(I(-4)))*rb)*(pow(ab,2)-pow(bd,2)-pow(ad,2)))*pow(rc,2)+((pow(bc,2)*(I(-8))+pow(ab,2)*(I(-8)))*pow(rd,3))*r+(pow(r,2)*(I(-4)))*(pow(bd,4)+(pow(ab,2)-pow(ad,2))*pow(bc,2)-pow(bd,2)*(pow(bc,2)+pow(ab,2)+pow(ad,2))-pow(cd,2)*(pow(bd,2)+pow(ab,2)-pow(ad,2)))+((pow(ab,2)*(I(-2)))*pow(bd,2))*pow(cd,2)+(pow(bc,2)*(I(-2))+pow(ab,2)*(I(-2)))*pow(rd,4)+(((pow(bd,2)*I(4)+pow(bc,2)*(I(-4))+pow(cd,2)*I(4))*pow(ra,2))*r+((pow(bd,2)*I(8)+pow(bc,2)*(I(-8))+pow(cd,2)*I(8))*ra)*pow(r,2)+((pow(bd,2)*I(4)+pow(ab,2)*I(4)+pow(ad,2)*(I(-8)))*pow(cd,2)+pow(bd,2)*I(4)*pow(ad,2)+pow(bc,2)*I(4)*pow(ad,2))*r)*rb+((pow(bd,2)*I(4)*pow(cd,2)+pow(bc,2)*I(4)*pow(bd,2)+pow(bd,4)*(I(-4)))*ra)*r+pow(rd,2)*I(2)*(((pow(bd,2)*I(2)+pow(ab,2)*I(2)+pow(ad,2)*(I(-2)))*r)*rc+ra*I(2)*r*(pow(bd,2)+pow(bc,2)-pow(cd,2))+pow(rb,2)*(pow(bc,2)+pow(ab,2)+pow(bd,2)*(I(-2))+pow(ad,2)+pow(cd,2))+((pow(bc,2)*I(2)+pow(ab,2)*I(2)+pow(bd,2)*(I(-4))+pow(ad,2)*I(2)+pow(cd,2)*I(2))*r)*rb+pow(cd,2)*pow(ab,2)+(pow(bc,2)*(I(-4))+pow(ab,2)*(I(-4)))*pow(r,2)+pow(rc,2)*(pow(bd,2)+pow(ab,2)-pow(ad,2))+(pow(bd,2)+pow(bc,2)-pow(cd,2))*pow(ra,2)+pow(bd,2)*(pow(bc,2)+pow(ab,2))-(pow(ab,2)*I(2)-pow(ad,2))*pow(bc,2))+(pow(bd,2)*I(2)*pow(cd,2)+pow(bc,2)*I(2)*pow(bd,2)+pow(bd,4)*(I(-2)))*pow(ra,2)+pow(rb,2)*I(2)*(((pow(bd,2)*I(2)+pow(bc,2)*(I(-2))+pow(cd,2)*I(2))*ra)*r+(pow(cd,2)*(I(-4))+pow(ad,2)*(I(-4)))*pow(r,2)+pow(ad,2)*pow(bd,2)+pow(cd,2)*(pow(bd,2)+pow(ab,2)+pow(ad,2)*(I(-2)))+pow(ad,2)*pow(bc,2)-pow(ra,2)*(pow(bc,2)-pow(bd,2)-pow(cd,2)))+((pow(bc,2)*(I(-2)))*pow(bd,2))*pow(ad,2);
    I F=r*I(4)*pow(rc,3)*(pow(ab,4)+(pow(ab,2)*(I(-2)))*pow(ad,2)+(pow(ad,2)*(I(-2))+pow(ab,2)*(I(-2)))*pow(bd,2)+pow(bd,4)+pow(ad,4))+pow(ra,3)*I(4)*r*(pow(bc,4)+(pow(bc,2)*(I(-2)))*pow(bd,2)+pow(bd,4)+(pow(bd,2)*(I(-2))+pow(bc,2)*(I(-2)))*pow(cd,2)+pow(cd,4))+(pow(ab,4)+(pow(ab,2)*(I(-2)))*pow(bc,2)+pow(bc,4))*pow(rd,4)+(ra*I(8)*pow(r,2)*(pow(ad,2)*pow(bc,2)-pow(ad,2)*pow(bd,2)-pow(cd,2)*(pow(ab,2)*I(2)-pow(bc,2)-pow(bd,2)-pow(ad,2))-pow(cd,4))+pow(ra,2)*I(4)*r*(pow(ad,2)*pow(bc,2)-pow(ad,2)*pow(bd,2)-pow(cd,2)*(pow(ab,2)*I(2)-pow(bc,2)-pow(bd,2)-pow(ad,2))-pow(cd,4))+(r*(I(-4)))*(pow(ad,4)*pow(bc,2)+pow(cd,4)*pow(ab,2)-pow(cd,2)*(pow(ad,2)*pow(bc,2)+pow(ad,2)*pow(ab,2))))*rb+((pow(bd,2)*(I(-4)))*(pow(ad,2)*pow(ab,2)-pow(ad,2)*pow(bc,2))+(pow(cd,2)*(I(-4)))*(pow(ab,4)-pow(ad,2)*pow(ab,2)-pow(bd,2)*(pow(ab,2)-pow(bc,2))-pow(bc,2)*(pow(ad,2)+pow(ab,2)))+(pow(ad,4)*(I(-4))+pow(ab,2)*I(4)*pow(ad,2))*pow(bc,2)+(pow(ab,2)*(I(-4)))*pow(cd,4)+(pow(bc,4)*(I(-4)))*pow(ad,2))*pow(r,2)+pow(rd,3)*I(4)*r*(pow(ab,4)+(pow(ab,2)*(I(-2)))*pow(bc,2)+pow(bc,4))+(((r*(I(-4)))*pow(rc,2))*(pow(ab,2)*I(2)*pow(cd,2)+pow(ab,4)-pow(ad,2)*pow(ab,2)-pow(bd,2)*(pow(ab,2)-pow(bc,2))-pow(bc,2)*(pow(ad,2)+pow(ab,2)))+((pow(r,2)*(I(-8)))*rc)*(pow(ab,2)*I(2)*pow(cd,2)+pow(ab,4)-pow(ad,2)*pow(ab,2)-pow(bd,2)*(pow(ab,2)-pow(bc,2))-pow(bc,2)*(pow(ad,2)+pow(ab,2)))+((ra*(I(-8)))*pow(r,2))*(pow(bc,4)+pow(bd,2)*(pow(ab,2)-pow(bc,2))-(pow(ad,2)*(I(-2))+pow(ab,2))*pow(bc,2)-pow(cd,2)*(pow(bc,2)+pow(ab,2)))+((pow(ra,2)*(I(-4)))*r)*(pow(bc,4)+pow(bd,2)*(pow(ab,2)-pow(bc,2))-(pow(ad,2)*(I(-2))+pow(ab,2))*pow(bc,2)-pow(cd,2)*(pow(bc,2)+pow(ab,2)))+((r*(I(-4)))*pow(rb,2))*(pow(ad,2)*pow(ab,2)-pow(cd,2)*(pow(ab,2)-pow(bc,2))-pow(ad,2)*pow(bc,2))+((pow(r,2)*(I(-8)))*rb)*(pow(ad,2)*pow(ab,2)-pow(cd,2)*(pow(ab,2)-pow(bc,2))-pow(ad,2)*pow(bc,2))+r*I(4)*(pow(bc,2)*pow(ab,2)*pow(ad,2)-(pow(ab,4)-pow(bc,2)*pow(ab,2))*pow(cd,2)-pow(ad,2)*pow(bc,4)))*rd+r*I(4)*pow(rb,3)*(pow(ad,4)+(pow(ad,2)*(I(-2)))*pow(cd,2)+pow(cd,4))+(pow(rb,2)*(I(-2)))*((pow(r,2)*(I(-2)))*(pow(ad,4)+(pow(ad,2)*(I(-2)))*pow(cd,2)+pow(cd,4))+((ra*(I(-2)))*r)*(pow(ad,2)*pow(bc,2)-pow(ad,2)*pow(bd,2)-pow(cd,2)*(pow(ab,2)*I(2)-pow(bc,2)-pow(bd,2)-pow(ad,2))-pow(cd,4))+pow(cd,4)*pow(ab,2)+pow(ad,4)*pow(bc,2)-pow(cd,2)*(pow(ad,2)*pow(bc,2)+pow(ad,2)*pow(ab,2))-(pow(ad,2)*pow(bc,2)-pow(ad,2)*pow(bd,2)-pow(cd,2)*(pow(ab,2)*I(2)-pow(bc,2)-pow(bd,2)-pow(ad,2))-pow(cd,4))*pow(ra,2))+(pow(ad,4)+(pow(ad,2)*(I(-2)))*pow(cd,2)+pow(cd,4))*pow(rb,4)+((pow(bd,2)*(I(-2)))*(pow(ab,2)*I(2)*pow(ad,2)-pow(ad,2)*pow(bc,2))+pow(rb,2)*I(2)*(pow(ad,2)*pow(bd,2)+pow(ad,2)*pow(ab,2)+pow(cd,2)*(pow(ad,2)+pow(ab,2)-pow(bd,2))+(pow(bc,2)*(I(-2)))*pow(ad,2)-pow(ad,4))+(pow(ad,4)*(I(-2))+pow(ab,2)*I(2)*pow(ad,2))*pow(bc,2)+(pow(cd,2)*(I(-2)))*(pow(ab,4)-pow(ad,2)*pow(ab,2)-pow(bd,2)*pow(ab,2))+(pow(ra,2)*(I(-2)))*(pow(bd,4)+(pow(ab,2)-pow(ad,2))*pow(bc,2)-pow(bd,2)*(pow(bc,2)+pow(ab,2)+pow(ad,2))-pow(cd,2)*(pow(bd,2)+pow(ab,2)-pow(ad,2)))+((ra*(I(-4)))*r)*(pow(bd,4)+(pow(ab,2)-pow(ad,2))*pow(bc,2)-pow(bd,2)*(pow(bc,2)+pow(ab,2)+pow(ad,2))-pow(cd,2)*(pow(bd,2)+pow(ab,2)-pow(ad,2)))+r*I(4)*rb*(pow(ad,2)*pow(bd,2)+pow(ad,2)*pow(ab,2)+pow(cd,2)*(pow(ad,2)+pow(ab,2)-pow(bd,2))+(pow(bc,2)*(I(-2)))*pow(ad,2)-pow(ad,4))+pow(r,2)*I(4)*(pow(ab,4)+(pow(ab,2)*(I(-2)))*pow(ad,2)+(pow(ad,2)*(I(-2))+pow(ab,2)*(I(-2)))*pow(bd,2)+pow(bd,4)+pow(ad,4)))*pow(rc,2)+(pow(ra,2)*(I(-2)))*(pow(cd,4)*pow(ab,2)+(pow(r,2)*(I(-2)))*(pow(bc,4)+(pow(bc,2)*(I(-2)))*pow(bd,2)+pow(bd,4)+(pow(bd,2)*(I(-2))+pow(bc,2)*(I(-2)))*pow(cd,2)+pow(cd,4))+pow(ad,2)*pow(bc,4)-pow(cd,2)*(pow(bc,2)*(pow(ad,2)+pow(ab,2))+pow(bd,2)*(pow(bc,2)*(I(-2))+pow(ab,2)))-pow(bd,2)*pow(bc,2)*pow(ad,2))+(((pow(ab,2)*(I(-2)))*pow(bc,2))*pow(ad,2))*pow(cd,2)+((ra*(I(-4)))*(pow(ad,2)*pow(bc,4)+pow(cd,4)*pow(ab,2)-pow(cd,2)*(pow(bc,2)*(pow(ad,2)+pow(ab,2))+pow(bd,2)*(pow(bc,2)*(I(-2))+pow(ab,2)))-pow(bd,2)*pow(bc,2)*pow(ad,2)))*r+pow(cd,4)*pow(ab,4)+(pow(bc,4)+(pow(bc,2)*(I(-2)))*pow(bd,2)+pow(bd,4)+(pow(bd,2)*(I(-2))+pow(bc,2)*(I(-2)))*pow(cd,2)+pow(cd,4))*pow(ra,4)+(pow(ab,4)+(pow(ab,2)*(I(-2)))*pow(ad,2)+(pow(ad,2)*(I(-2))+pow(ab,2)*(I(-2)))*pow(bd,2)+pow(bd,4)+pow(ad,4))*pow(rc,4)+pow(rd,2)*I(2)*(pow(r,2)*I(2)*(pow(ab,4)+(pow(ab,2)*(I(-2)))*pow(bc,2)+pow(bc,4))+((r*(I(-2)))*rc)*(pow(ab,2)*I(2)*pow(cd,2)+pow(ab,4)-pow(ad,2)*pow(ab,2)-pow(bd,2)*(pow(ab,2)-pow(bc,2))-pow(bc,2)*(pow(ad,2)+pow(ab,2)))+((ra*(I(-2)))*r)*(pow(bc,4)+pow(bd,2)*(pow(ab,2)-pow(bc,2))-(pow(ad,2)*(I(-2))+pow(ab,2))*pow(bc,2)-pow(cd,2)*(pow(bc,2)+pow(ab,2)))+((r*(I(-2)))*rb)*(pow(ad,2)*pow(ab,2)-pow(cd,2)*(pow(ab,2)-pow(bc,2))-pow(ad,2)*pow(bc,2))+pow(bc,2)*pow(ab,2)*pow(ad,2)-(pow(bc,4)+pow(bd,2)*(pow(ab,2)-pow(bc,2))-(pow(ad,2)*(I(-2))+pow(ab,2))*pow(bc,2)-pow(cd,2)*(pow(bc,2)+pow(ab,2)))*pow(ra,2)-(pow(ab,4)-pow(bc,2)*pow(ab,2))*pow(cd,2)-(pow(ad,2)*pow(ab,2)-pow(cd,2)*(pow(ab,2)-pow(bc,2))-pow(ad,2)*pow(bc,2))*pow(rb,2)-(pow(ab,2)*I(2)*pow(cd,2)+pow(ab,4)-pow(ad,2)*pow(ab,2)-pow(bd,2)*(pow(ab,2)-pow(bc,2))-pow(bc,2)*(pow(ad,2)+pow(ab,2)))*pow(rc,2)-pow(ad,2)*pow(bc,4))+pow(ad,4)*pow(bc,4)+(rc*(I(-4)))*(r*pow(ra,2)*(pow(bd,4)+(pow(ab,2)-pow(ad,2))*pow(bc,2)-pow(bd,2)*(pow(bc,2)+pow(ab,2)+pow(ad,2))-pow(cd,2)*(pow(bd,2)+pow(ab,2)-pow(ad,2)))+ra*I(2)*pow(r,2)*(pow(bd,4)+(pow(ab,2)-pow(ad,2))*pow(bc,2)-pow(bd,2)*(pow(bc,2)+pow(ab,2)+pow(ad,2))-pow(cd,2)*(pow(bd,2)+pow(ab,2)-pow(ad,2)))+((pow(r,2)*(I(-2)))*rb)*(pow(ad,2)*pow(bd,2)+pow(ad,2)*pow(ab,2)+pow(cd,2)*(pow(ad,2)+pow(ab,2)-pow(bd,2))+(pow(bc,2)*(I(-2)))*pow(ad,2)-pow(ad,4))-(pow(bc,2)*(pow(ad,2)*pow(ab,2)-pow(ad,4))-(pow(ab,2)*I(2)*pow(ad,2)-pow(ad,2)*pow(bc,2))*pow(bd,2)-(pow(ab,4)-pow(ad,2)*pow(ab,2)-pow(bd,2)*pow(ab,2))*pow(cd,2))*r-pow(rb,2)*r*(pow(ad,2)*pow(bd,2)+pow(ad,2)*pow(ab,2)+pow(cd,2)*(pow(ad,2)+pow(ab,2)-pow(bd,2))+(pow(bc,2)*(I(-2)))*pow(ad,2)-pow(ad,4)));




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
    if (!zero_in(D*pow(H,4)+E*pow(H,2)+F)) {printf("no way\n");return I(1,0);}

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
        if (upper(x)<0.78) acs2=intersect(acs2,-F/E*(I(1)+x));
    }


    I acs=sqrt(acs2);
    acs=intersect(acs,H);
    return acs;
}

