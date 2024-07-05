// the radius is the root of A*x^2+B*x*C where A,B,C are complicate polynomials in the edge lengths and sphere radii
// the formulas for A,B,C are here specialized, depending on the radii and assuming that there is contact along ab & ac
// cases given by the radii a,b,c: uuu, uur, urr, ruu, rur, rrr
// TODO : simplified formulas for these cases!

I radius(I ab, I ac, I ad, I bc, I bd, I cd, bool verbose=false)
{
    I w=square(ad);
    I x=square(bc);
    I y=square(bd);
    I z=square(cd);
    #if defined (uuuu)
    I A=I(4)*x*(w*(I(8)+y+z-x-w)+(z-I(4))*(I(4)-y))-I(16)*square(y-z);
    I B=I(2)*A;
    I C=x*((w-I(4))*(w*(x-I(4)) - I(4)*y - I(4)*z+I(16)) - I(4)*y*z);
#else//if defined (ruuu)
    I A=(square(I(-2)+w))*((I(-4))*x)+(w+I(-2)+I(8)*r)*(y*I(4)*x)+r*((I(-16))*square(y))+r*((I(-16))*square(z))+(w+I(-2)+I(4)*r)*((I(-4))*square(x))+(x*(w+I(-2)+I(8)*r)+y*(I(8)*r-x))*(I(4)*z);
    I B=(x*((I(1)+r)*(w)+I(6)+(I(-10))*r)+(I(4)*r+I(-4)+x*r)*((I(-2))*y))*(I(4)*z)+(I(-1)+r)*(I(16)*square(y))+(square(I(-2)+w))*((I(-8))*x)+(I(-1)+r)*(I(16)*square(z))+((I(1)+r)*(w)+I(6)+(I(-10))*r)*(y*I(4)*x)+((I(1)+r)*(w)+I(2)+(I(-6))*r)*((I(-4))*square(x));
    I C=square(x)*(square(w)+(I(-8))*r+I(4)+(I(-1)+r)*(I(4)*w))+(square(I(-2)+w))*((I(-4))*x)+((I(-2))*r+r*w)*(y*((I(-4))*x))+((I(-1)+I(2)*r)*(x*y)-x*((I(-2))*r+r*w))*(I(4)*z);
    #endif
    if (verbose)
    {
        printf("A=%f,%f\n",lower(A),upper(A));
        printf("B=%f,%f\n",lower(B),upper(B));
        printf("C=%f,%f\n",lower(C),upper(C));
        I X=(A*hull(I(0),r)+B)*hull(I(0),r)+C;
        printf("A*r^2+B*r+C=%f,%f\n",lower(X),upper(X));
    }
    // no chance to have a root smaller than r
    if (!zero_in((A*hull(I(0),r)+B)*hull(I(0),r)+C)) return I(INFINITY);
    // usual formula through discriminant
    I D=square(B)-I(4)*A*C;
    if (upper(D)<=0) return I(INFINITY); // no real root (seems never happen)
    I rs=(upper(C)<0 ? (-B+sqrt(D))/(I(2)*A) : (-B-sqrt(D))/(I(2)*A));
    // if A contains 0 or is very small: use the continuity of the roots in the coefficients!
    if (zero_in(A)||zero_in(C))
    {
        I X=I(4)*A*C/square(B);
        X=hull(x,I(0)); // the interval x must contain 0
        if (upper(x)<0.78) rs=intersect(rs,-C/B*(I(1)+X));
    }
    return rs;
}

