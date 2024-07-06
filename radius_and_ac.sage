
# This file shows how to find an explicit formula for the radius of the support sphere of a tetrahedron

# known variables: edge length and sphere radii
var('ab ac ad bc bd cd ra rb rc rd')

# first, we express the coordinates of the vertices of the tetrahedron in the known variables
# we introduce new unknown variables:
var('xa ya za xb yb zb xc yc zc xd yd zd')
# the dictionnary val will collect expressions in the known variables for these new variables
val={}
# we set the axis so that A=(0,0,0), B=(ab,0,0), C=(xc,yc,0), yc>0 and zd>0
val.update({xa:0,ya:0,za:0,xb:ab,yb:0,zb:0,zc:0})
# the following equations relate vertex coordinates and edge length
eq_ab=(xb-xa)^2+(yb-ya)^2+(zb-za)^2==ab^2
eq_ac=(xc-xa)^2+(yc-ya)^2+(zc-za)^2==ac^2
eq_ad=(xd-xa)^2+(yd-ya)^2+(zd-za)^2==ad^2
eq_bc=(xc-xb)^2+(yc-yb)^2+(zc-zb)^2==bc^2
eq_bd=(xd-xb)^2+(yd-yb)^2+(zd-zb)^2==bd^2
eq_cd=(xd-xc)^2+(yd-yc)^2+(zd-zc)^2==cd^2
# let us find an expression
# xc can be computed with eq_ac and eq_bc:
val.update({xc:solve(eq_bc-eq_ac,xc)[0].right_hand_side().subs(val)})
# yc is then the positive solution of eq_ac
val.update({yc:solve(eq_ac,yc)[1].right_hand_side().subs(val)})
# xd can then be computed with eq_ad and eq_bd:
val.update({xd:solve(eq_bd-eq_ad,xd)[0].right_hand_side().subs(val)})
# yd can then be computed with eq_bd and eq_cd:
val.update({yd:solve(eq_bd-eq_cd,yd)[0].right_hand_side().subs(val)})
# zd can then be computed with eq_ad
val.update({zd:solve(eq_ad,zd)[0].right_hand_side().subs(val)})


# now, we express the center (x,y,z) and radius R of the support circle in the known variables
var('x y z R')
# the following equations follow from the definition of the support circle
eq1=(xa-x)^2+(ya-y)^2+(za-z)^2==(R+ra)^2
eq2=(xb-x)^2+(yb-y)^2+(zb-z)^2==(R+rb)^2
eq3=(xc-x)^2+(yc-y)^2+(zc-z)^2==(R+rc)^2
eq4=(xd-x)^2+(yd-y)^2+(zd-z)^2==(R+rd)^2
# we can easily get three LINEAR equations in x,y,z,R:
eq5=expand(eq1-eq2)
eq6=expand(eq1-eq3)
eq7=expand(eq1-eq4)
# this allows to express x,y,z with R and the coordinates of the vertices:
t=solve([eq5,eq6,eq7],[x,y,z])[0]
# the dictionnary val2 collects expressions for x,y,z in the known variables
val2={
x:t[0].right_hand_side().subs(val).full_simplify(),
y:t[1].right_hand_side().subs(val).full_simplify(),
z:t[2].right_hand_side().subs(val).full_simplify()
}
# if we replace x,y,z by their expression in eq1, we get a quadratic polynomial in R and the known variables.

# Denote this polynomial by P=A*R^2+B*R+C and let us find an expression for its coeffients:
P=(xa-x)^2+(ya-y)^2+(za-z)^2-(R+ra)^2
P=P.subs(val).subs(val2).expand()
P=numerator(P)
A=P.coefficient(R^2).expand()
B=P.coefficient(R).expand()
C=P.subs(R=0).expand()
# the coefficients A,B,C are "not very nice" but are polynomials in the known variables

# R is the root (-B-sign(C)*sqrt(D))/(2A)

var('r')

# formula for the length of ac when the support sphere is assumed to have radius r
P=PolynomialRing(QQ[ab,ad,bc,bd,cd,ra,rb,rc,rd,r],ac)(A*r^2+B*r+C)
# P is a quadratic polynomials in ac^2
D=SR(P).coefficient(ac^4).expand()
E=SR(P).coefficient(ac^2).expand()
F=SR(P).subs(ac=0).expand()

# D=(bd + 2*r + rb + rd)*(bd - 2*r - rb - rd)*(bd + rb - rd)*(bd - rb + rd)
# D is always nonpositive

# ac^2 is equal to -E/2/D plus or minus sqrt(E^2-4*D*F)/2/D
# if the discriminant is negative or there is no solution in [ra+rb,ra+rb+ab], then there is no length ac which allows such a FM-tetrahedron
# what if two valid solutions? may happen?

# full_simplify E and F, replace ab^2 by ab2
var('u x y z w')
dico_ro = {ad^2:w, ac^2:u, cd^2:z, bd^2:y, bc^2:x, ad^4:w^2, ac^4:u^2, cd^4:z^2, bd^4:y^2, bc^4:x^2}
dico_1111_ab_cd = {ra:1,rb:1,rc:1,rd:1, ab:2, cd:2}
rv = sqrt(2)-1
dico_rrrr_ab_cd = {ra:r,rb:r,rc:r,rd:r, ab:2*r, cd:2*r}
dico_rrrr_ab_ac = {ra:r,rb:r,rc:r,rd:r, ab:2*r, ac:2*r}
dico_rrrr_ab = {ra:r,rb:r,rc:r,rd:r, ab:2*r}

AB,AC,AD,BC,BD,CD = (RIF(0.82842712474618984686003386741504073143005371093750,0.82842712474619029094924371747765690088272094726562),
                     RIF(0.82842712474618984686003386741504073143005371093750,1.65685424949238058189848743495531380176544189453125),
                     RIF(0.82842713594436645507812500000000000000000000000000,0.83878247439861297607421875000000000000000000000000),
                     RIF(1.09766596555709838867187500000000000000000000000000,1.11837661266326904296875000000000000000000000000000),
                     RIF(1.59472227096557617187500000000000000000000000000000,1.60507762432098388671875000000000000000000000000000),
		     RIF(0.82842713594436645507812500000000000000000000000000,0.84913781285285949707031250000000000000000000000000) )
