// the radius is the root of A*x^2+B*x*C where A,B,C are complicate polynomials in the edge lengths and sphere radii
// the formulas for A,B,C are here specialized, depending on the radii and assuming some contacts
// cases given by the radii a,b,c: uuu, uur, urr, ruu, rur, rrr
// TODO : simplified formulas for these cases!

I radius(I ab, I ac, I ad, I bc, I bd, I cd, bool verbose=false)
{
    I x=square(bc);
    I y = square(bd);
    I z = square(cd);
    I w = square(ad);
    I u = square(ac);
#if defined (uuuu) && defined(contact_ac)
    I A=I(4)*x*(w*(I(8)+y+z-x-w)+(z-I(4))*(I(4)-y))-I(16)*square(y-z);
    I B=I(2)*A;
    I C=x*((w-I(4))*(w*(x-I(4)) - I(4)*y - I(4)*z+I(16)) - I(4)*y*z);
#elif defined(uuuu) && defined(contact_cd)
    I A = I(4)*((u*(I(-4)))*x+I(-128)+y*I(16)+(y*(I(-4)))*w+w*I(8)*x+y*I(8)*u+w*u*x+u*y*w+u*I(16)+w*I(16)+x*I(16)+u*y*x+(u*(I(-4)))*w+(y*(I(-4)))*x+w*y*x-x*pow(w,2)-pow(x,2)*w-u*pow(y,2)-pow(u,2)*y);
    I B = I(2)*A;
    I C = (u*(I(-16)))*x+I(-256)+y*I(64)+(y*(I(-16)))*w+(pow(w,2)*(I(-4)))*x+(w*(I(-4)))*pow(x,2)+pow(x,2)*pow(w,2)+(pow(y,2)*(I(-4)))*u+pow(u,2)*pow(y,2)+(y*(I(-4)))*pow(u,2)+u*I(4)*w*x+y*I(4)*u*w+u*I(64)+w*I(64)+(((y*(I(-2)))*u)*w)*x+x*I(64)+y*I(4)*u*x+(u*(I(-16)))*w+(y*(I(-16)))*x+y*I(4)*w*x;
#elif defined(rrrr) && defined(contact_ac) // TODO check
    I A = I(4)*(square(r)*z*I(4)*x+square(r)*y*I(8)*z+square(r)*x*I(8)*w+y*w*x+square(r)*((I(-4))*square(z))+z*w*x+square(r)*y*I(4)*x+pow(r,4)*((I(-16))*x)+square(r)*((I(-4))*square(y))-square(w)*x-w*square(x)-y*x*z);
    I B = I(4)*(pow(r,3)*z*I(8)*x+pow(r,3)*y*I(16)*z+pow(r,3)*x*I(16)*w+r*(square(x)*((I(-2))*w))+r*(x*((I(-2))*square(w)))+pow(r,3)*((I(-8))*square(z))+pow(r,3)*y*I(8)*x+pow(r,3)*((I(-8))*square(y))+r*(y*(z*((I(-2))*x)))+r*z*x*I(2)*w+r*y*x*I(2)*w+pow(r,5)*((I(-32))*x));
    I C = pow(r,4)*z*I(16)*x+square(r)*(x*((I(-4))*square(w)))+square(r)*(square(x)*((I(-4))*w))+pow(r,4)*x*I(32)*w+pow(r,4)*y*I(16)*x+square(w)*square(x)+square(r)*(y*(z*((I(-4))*x)))+square(r)*(z*(x*((I(-4))*w)))+square(r)*(y*(x*((I(-4))*w)))+pow(r,6)*((I(-64))*x) ;
#elif defined(rrrr) && defined(contact_cd) // TODO check
    I A = I(4)* ((I(-128))*pow(r,6)+u*x*y+square(r)*x*I(8)*w+y*w*x+u*w*y+square(r)*(y*((I(-4))*x))+u*I(16)*pow(r,4)+u*w*x+square(r)*(y*((I(-4))*w))+pow(r,4)*I(16)*y+u*square(r)*I(8)*y+u*(square(r)*((I(-4))*x))+pow(r,4)*I(16)*x+pow(r,4)*I(16)*w+u*(square(r)*((I(-4))*w))-square(w)*x-w*square(x)-square(y)*u-y*square(u));
    I B = I(4)*((I(-256))*pow(r,7)+pow(r,3)*x*I(16)*w+r*(square(x)*((I(-2))*w))+r*(x*((I(-2))*square(w)))+u*r*x*I(2)*w+pow(r,3)*(y*((I(-8))*x))+u*I(32)*pow(r,5)+u*r*y*I(2)*w+pow(r,3)*(y*((I(-8))*w))+pow(r,5)*I(32)*y+u*(r*((I(-2))*square(y)))+u*pow(r,3)*I(16)*y+square(u)*(r*((I(-2))*y))+u*r*y*I(2)*x+u*(pow(r,3)*((I(-8))*x))+pow(r,5)*I(32)*x+u*(pow(r,3)*((I(-8))*w))+pow(r,5)*I(32)*w+r*y*x*I(2)*w);
    I C =(I(-256))*pow(r,8)+square(y)*square(u)+square(r)*(x*((I(-4))*square(w)))+square(r)*(square(x)*((I(-4))*w))+u*(y*(x*((I(-2))*w)))+u*square(r)*x*I(4)*w+pow(r,4)*(y*((I(-16))*x))+u*I(64)*pow(r,6)+square(w)*square(x)+u*square(r)*y*I(4)*w+pow(r,4)*(y*((I(-16))*w))+pow(r,6)*I(64)*y+u*(square(r)*((I(-4))*square(y)))+square(u)*(square(r)*((I(-4))*y))+u*square(r)*y*I(4)*x+u*(pow(r,4)*((I(-16))*x))+pow(r,6)*I(64)*x+u*(pow(r,4)*((I(-16))*w))+pow(r,6)*I(64)*w+square(r)*y*x*I(4)*w;
#else//TODO other cases
    printf("Error: radius formula are not implemented for this case");
	I A = I(0);
    I B = I(0);
    I C = I(0);
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

