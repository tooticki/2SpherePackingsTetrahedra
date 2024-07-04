

# radius of the support sphere (if any) of the input tetrahedron
def radius(ab,ac,ad,bc,bd,cd,ra,rb,rc,rd,dessine=false):
    # coordinates of vertices
    # A:(0,0,0), B:(ab,0,0), C:(xc,yc>0,0), D(xd,yd,zd>0)
    rr = sqrt(2)-1
    xc=(ab^2+ac^2-bc^2)/(2*ab)
    yc=ac^2-xc^2
    if parent(yc)==RIF:
        yc=yc.intersection(RIF(0,infinity)) # si en  RIF
    yc=sqrt(yc)
    xd=(ab^2+ad^2-bd^2)/(2*ab)
    yd=RIF(-1/2)*(ab^2 - bd^2 + cd^2 - xc^2 - 2*ab*xd + 2*xc*xd - yc^2)/yc
    zd=ad^2-xd^2-yd^2
    if parent(zd)==RIF:
        zd=zd.intersection(RIF(0,infinity)) # si en  RIF
    zd=sqrt(zd)
    # draw the tetrahedron

    # radius of support sphere
    A=4*(ab^2*ra^2 - 2*ab^2*ra*rd + ab^2*rd^2 + (ra^2 - 2*ra*rb + rb^2)*xd^2 - 2*(ab*ra^2 - ab*ra*rb - (ab*ra - ab*rb)*rd)*xd)*yc^2 - 8*(ab^2*ra^2 - ab^2*ra*rc - (ab^2*ra - ab^2*rc)*rd - (ab*ra^2 - ab*ra*rb - (ab*ra - ab*rb)*rd)*xc - (ab*ra^2 - ab*ra*rb - (ab*ra - ab*rb)*rc - (ra^2 - 2*ra*rb + rb^2)*xc)*xd)*yc*yd + 4*(ab^2*ra^2 - 2*ab^2*ra*rc + ab^2*rc^2 + (ra^2 - 2*ra*rb + rb^2)*xc^2 - 2*(ab*ra^2 - ab*ra*rb - (ab*ra - ab*rb)*rc)*xc)*yd^2 + 4*(ab^2*ra^2 - 2*ab^2*ra*rc + ab^2*rc^2 + (ra^2 - 2*ra*rb + rb^2)*xc^2 - (ab^2 - ra^2 + 2*ra*rb - rb^2)*yc^2 - 2*(ab*ra^2 - ab*ra*rb - (ab*ra - ab*rb)*rc)*xc)*zd^2
    B=-4*(ab^2*ra - ab^2*rc - (ab*ra - ab*rb)*xc)*yc*yd^3 + 4*(ab^2*ra^3 - ab^2*ra^2*rd - ab^2*ra*rd^2 + ab^2*rd^3 - (ab*ra - ab*rb)*xd^3 + (2*ab^2*ra + ra^3 - ra*rb^2 + rb^3 - ab^2*rd - (ab^2 + ra^2)*rb)*xd^2 - (ab^3*ra + 2*ab*ra^3 - ab*ra^2*rb - ab*ra*rb^2 - (ab*ra - ab*rb)*rd^2 - (ab^3 + ab*ra^2 - ab*rb^2)*rd)*xd)*yc^2 + 4*(ab^2*ra^3 - ab^2*ra^2*rc - ab^2*ra*rc^2 + ab^2*rc^3 - (ab*ra - ab*rb)*xc^3 + (2*ab^2*ra + ra^3 - ra*rb^2 + rb^3 - ab^2*rc - (ab^2 + ra^2)*rb)*xc^2 + (2*ab^2*ra - ab^2*rc - ab^2*rd - (ab*ra - ab*rb)*xc - (ab*ra - ab*rb)*xd)*yc^2 - (ab^3*ra + 2*ab*ra^3 - ab*ra^2*rb - ab*ra*rb^2 - (ab*ra - ab*rb)*rc^2 - (ab^3 + ab*ra^2 - ab*rb^2)*rc)*xc)*yd^2 + 4*(ab^2*ra^3 - ab^2*ra^2*rc - ab^2*ra*rc^2 + ab^2*rc^3 - (ab*ra - ab*rb)*xc^3 + (2*ab^2*ra + ra^3 - ra*rb^2 + rb^3 - ab^2*rc - (ab^2 + ra^2)*rb)*xc^2 + (ab^2*ra + ra^3 - ra*rb^2 + rb^3 - ab^2*rc - ab^2*rd - (ab^2 + ra^2)*rb - (ab*ra - ab*rb)*xc - (ab*ra - ab*rb)*xd)*yc^2 - (ab^2*ra - ab^2*rc - (ab*ra - ab*rb)*xc)*yc*yd - (ab^3*ra + 2*ab*ra^3 - ab*ra^2*rb - ab*ra*rb^2 - (ab*ra - ab*rb)*rc^2 - (ab^3 + ab*ra^2 - ab*rb^2)*rc)*xc)*zd^2 - 4*((ab^2*ra - ab^2*rd - (ab*ra - ab*rb)*xd)*yc^3 + (2*ab^2*ra^3 - ab^2*ra^2*rc - ab^2*ra*rc^2 - (ab^2*ra - ab^2*rc)*rd^2 + (ab^2*ra - ab^2*rd)*xc^2 + (ab^2*ra - ab^2*rc - (ab*ra - ab*rb)*xc)*xd^2 - (ab^2*ra^2 - ab^2*rc^2)*rd - (ab^3*ra + 2*ab*ra^3 - ab*ra^2*rb - ab*ra*rb^2 - (ab*ra - ab*rb)*rd^2 - (ab^3 + ab*ra^2 - ab*rb^2)*rd)*xc - (ab^3*ra + 2*ab*ra^3 - ab*ra^2*rb - ab*ra*rb^2 - (ab*ra - ab*rb)*rc^2 + (ab*ra - ab*rb)*xc^2 - (ab^3 + ab*ra^2 - ab*rb^2)*rc - 2*(ab^2*ra + ra^3 - ra*rb^2 + rb^3 - (ab^2 + ra^2)*rb)*xc)*xd)*yc)*yd
    C=ab^2*yc^2*yd^4 + ab^2*yc^2*zd^4 - 2*(ab^2*yc^3 + (ab^2*ra^2 - ab^2*rc^2 + ab^2*xc^2 - (ab^3 + ab*ra^2 - ab*rb^2)*xc)*yc)*yd^3 + (ab^2*ra^4 - 2*ab^2*ra^2*rd^2 + ab^2*rd^4 + ab^2*xd^4 - 2*(ab^3 + ab*ra^2 - ab*rb^2)*xd^3 + (ab^4 + 4*ab^2*ra^2 + ra^4 + rb^4 - 2*ab^2*rd^2 - 2*(ab^2 + ra^2)*rb^2)*xd^2 - 2*(ab^3*ra^2 + ab*ra^4 - ab*ra^2*rb^2 - (ab^3 + ab*ra^2 - ab*rb^2)*rd^2)*xd)*yc^2 + (ab^2*ra^4 - 2*ab^2*ra^2*rc^2 + ab^2*rc^4 + ab^2*xc^4 + ab^2*yc^4 - 2*(ab^3 + ab*ra^2 - ab*rb^2)*xc^3 + (ab^4 + 4*ab^2*ra^2 + ra^4 + rb^4 - 2*ab^2*rc^2 - 2*(ab^2 + ra^2)*rb^2)*xc^2 + 2*(2*ab^2*ra^2 - ab^2*rc^2 - ab^2*rd^2 + ab^2*xc^2 + ab^2*xd^2 - (ab^3 + ab*ra^2 - ab*rb^2)*xc - (ab^3 + ab*ra^2 - ab*rb^2)*xd)*yc^2 - 2*(ab^3*ra^2 + ab*ra^4 - ab*ra^2*rb^2 - (ab^3 + ab*ra^2 - ab*rb^2)*rc^2)*xc)*yd^2 + (ab^2*ra^4 - 2*ab^2*ra^2*rc^2 + ab^2*rc^4 + ab^2*xc^4 + ab^2*yc^4 + 2*ab^2*yc^2*yd^2 - 2*(ab^3 + ab*ra^2 - ab*rb^2)*xc^3 + (ab^4 + 4*ab^2*ra^2 + ra^4 + rb^4 - 2*ab^2*rc^2 - 2*(ab^2 + ra^2)*rb^2)*xc^2 + (ab^4 + 2*ab^2*ra^2 + ra^4 + rb^4 - 2*ab^2*rc^2 - 2*ab^2*rd^2 + 2*ab^2*xc^2 + 2*ab^2*xd^2 - 2*(ab^2 + ra^2)*rb^2 - 2*(ab^3 + ab*ra^2 - ab*rb^2)*xc - 2*(ab^3 + ab*ra^2 - ab*rb^2)*xd)*yc^2 - 2*(ab^3*ra^2 + ab*ra^4 - ab*ra^2*rb^2 - (ab^3 + ab*ra^2 - ab*rb^2)*rc^2)*xc - 2*(ab^2*yc^3 + (ab^2*ra^2 - ab^2*rc^2 + ab^2*xc^2 - (ab^3 + ab*ra^2 - ab*rb^2)*xc)*yc)*yd)*zd^2 - 2*((ab^2*ra^2 - ab^2*rd^2 + ab^2*xd^2 - (ab^3 + ab*ra^2 - ab*rb^2)*xd)*yc^3 + (ab^2*ra^4 - ab^2*ra^2*rc^2 - (ab^2*ra^2 - ab^2*rc^2)*rd^2 + (ab^2*ra^2 - ab^2*rd^2)*xc^2 + (ab^2*ra^2 - ab^2*rc^2 + ab^2*xc^2 - (ab^3 + ab*ra^2 - ab*rb^2)*xc)*xd^2 - (ab^3*ra^2 + ab*ra^4 - ab*ra^2*rb^2 - (ab^3 + ab*ra^2 - ab*rb^2)*rd^2)*xc - (ab^3*ra^2 + ab*ra^4 - ab*ra^2*rb^2 - (ab^3 + ab*ra^2 - ab*rb^2)*rc^2 + (ab^3 + ab*ra^2 - ab*rb^2)*xc^2 - (ab^4 + 2*ab^2*ra^2 + ra^4 + rb^4 - 2*(ab^2 + ra^2)*rb^2)*xc)*xd)*yc)*yd
    A=RIF(A)
    B=RIF(B)
    C=RIF(C)
    # no chance for a radius less than r to be root of Ax^2+Bx+C ! useful if A and B zero
    if not (A*RIF(0,rr)^2+B*RIF(0,rr)+C).contains_zero():
        return RIF(Infinity)
    # si r>=1 c'est vrai ; si r<=1: r^2<1 donc A*r^2+B*r+C=0 implique r>=-C/(A+B) et c'est donc encore vrai si -C/(A+B)>=1 ; utile
    #    if (-C/(A+B)).lower()>=1:
    #        return RIF(1,Infinity)
    D2=B^2-4*A*C # square of the discriminant
    if D2.upper()<=0: # no support sphere
        R=r+1 # any value larger than r will be recognized as invalid in a FM-tetrahedrization
    elif (A.contains_zero() or C.contains_zero()):
        x=RIF(0,RIF(4)*A*C/B^2)
        R=(-C/B*(1+x) if x.upper()<0.78 else RIF(-Infinity,Infinity))
    else:
        if parent(D2)==RIF:
            D2=D2.intersection(RIF(0,infinity))
            D=sqrt(D2)
            R=((-B+D)/(RIF(2)*A) if C<0 else (-B-D)/(RIF(2)*A))
    return R 

def dist2(x1,y1,z1,x2,y2,z2):
    return (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2

var('ab ad bc bd cd ra rb rc rd r')

#x1 == a*x+b*y+c, y1 == b*x-a*y+d
# 1->3 2->4
def find_transformation(x1,y1,x2,y2,x3,y3,x4,y4):
    a =((x1 - x2)*x3 - x1*x4 + x2*x4 - y1*(y3 - y4) + y2*y3 - y2*y4)/(x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2)
    b = ((x3 - x4)*y1 - x3*y2 + x4*y2 + x1*(y3 - y4) - x2*y3 + x2*y4)/(x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2)
    c = (x1^2*x4 + x4*y1^2 - (x2*x4 + y2*y3 - y2*y4)*x1 - (x1*x2 - x2^2 - y2^2)*x3 - (x3*y2 + x4*y2 - x2*y3 + x2*y4)*y1)/(x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2)
    d = (x1*x3*y2 + x1^2*y4 + y1^2*y4 - (x4*y2 + x2*y3 + x2*y4)*x1 - (x2*x3 - x2*x4 + y2*y3 + y2*y4)*y1 + (x2^2 + y2^2)*y3)/(x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2)
    return((a,b,c,d))
    

def point_by_three_distances(x0,y0,z0, x1,y1,z1,x2,y2,z2,d0,d1,d2):
    d01_2 = dist2(x0,y0,z0,x1,y1,z1)
    d02_2 = dist2(x0,y0,z0,x2,y2,z2)
    d12_2 = dist2(x2,y2,z2,x1,y1,z1)
    d01 = sqrt(d01_2)
    d02 = sqrt(d02_2)
    d12 = sqrt(d12_2)
    xx2 = (d02_2- d12_2 + d01_2)/2/d01
    yy2 = sqrt(d02_2 - xx2^2)
    print((xx2,yy2))
    x = (d0^2-d1^2+d01_2)/2/d01
    y = (d0^2-d2^2-2*xx2*x+xx2^2+yy2^2)/2/yy2
    z = sqrt(d0^2-x^2-y^2)
    # inverse to 0->000 1->0x0
    (a,b,c,d) = find_transformation(0,0,0,x0,y0,z0,0,d01,0,x1,y1,z1)
    return((x,y,z))




xo, yo, zo = 0,0,0
xa, ya, za = r+ra,0,0
xb = (-ab^2+(r+rb)^2+(r+ra)^2)/2/(r+ra)
yb = sqrt((r+rb^2)-xb^2)
zb = 0
xd,yd,zd = point_by_three_distances(xo,yo,zo,xa,ya,za,xb,yb,zb,r+rd,ad,bd)
xc,yc,zc = point_by_three_distances(xo,yo,zo,xd,yd,zd,xb,yb,zb,r+rc,cd,bc)
ac = sqrt((xc-r-ra)^2 + yc^2 + zc^2)

formula_ac = sqrt((2*r^4 + 2*r^3*ra + r^2*ra^2 + 2*r^3*rc + r^2*rc^2 + (2*r^2 + 2*r*ra + ra^2 + 2*r*rc + rc^2)*rd^2 + 2*(2*r^3 + 2*r^2*ra + r*ra^2 + 2*r^2*rc + r*rc^2)*rd + (cd^2*r - 2*r^3 - (r + ra)*rc^2 - (r + ra)*rd^2 + (cd^2 - 2*r^2)*ra - 2*(r^2 + r*ra)*rc - 2*(r^2 + r*ra)*rd)*sqrt(r^2 + 2*r*rd + rd^2))/(r^2 + 2*r*rd + rd^2))

def check (ab1,ac1,ad1,bc1,bd1,cd1,ra1,rb1,rc1,rd1):
    rr = radius(ab1,ac1,ad1,bc1,bd1,cd1,ra1,rb1,rc1,rd1)
    print (formula_ac.subs({ab:ab1,ad:ad1,bc:bc1,bd:bd1,cd:cd1,ra:ra1,rb:rb1,rc:rc1,rd:rd1,r:rr}))
    print(ac1)
    
u = RIF(1)
