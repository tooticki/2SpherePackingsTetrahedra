def sage2C(exp):
    (x,y)=(SR.wild(0),SR.wild(1))
    #print(exp)
    #print("  --  ")
    if exp.is_integer():
        return 'I('+str(exp)+')'
    if exp.is_constant():
        return 'I('+str(n(exp))+')'
    if exp.is_symbol():
        return str(exp)    

    m=exp.match(1/sqrt(y))
    if str(m)!='None':
        Y=sage2C(m[y])   
        return ('I(1) / sqrt('+Y+')')

    
    m=exp.match(x/y)
    if str(m)!='None':
        X=sage2C(m[x])
        if '+' in X or '-' in X: # ça serait mieux de déterminer si c'est vraiment un produit...
            X='('+str(X)+')'
        Y=sage2C(m[y])
        if '+' in Y or '-' in Y or '*' in Y:
            Y='('+str(Y)+')'    
        return (X+' / '+Y)

    m=exp.match(x*y)
    if str(m)!='None':
        X=sage2C(m[x])
        if '+' in X or '-' in X: # ça serait mieux de déterminer si c'est vraiment un produit...
            X='('+str(X)+')'
        Y=sage2C(m[y])
        if '+' in Y or '-' in Y:
            Y='('+str(Y)+')'
        return X+'*'+Y

    m=exp.match(x-y) # inutile ? il semble toujours matcher x-y en x+(-y)
    if str(m)!='None':
        return str(sage2C(m[x]))+'-'+str(sage2C(m[y]))
    m=exp.match(x+y)
    if str(m)!='None':
        return str(sage2C(m[x]))+'+'+str(sage2C(m[y]))
    m=exp.match(x^2)
    if str(m)!='None':
        return 'pow('+str(sage2C(m[x]))+',2)'
    m=exp.match(x^3)
    if str(m)!='None':
        return 'pow('+str(sage2C(m[x]))+',3)'
    m=exp.match(x^4)
    if str(m)!='None':
        return 'pow('+str(sage2C(m[x]))+',4)'
    m=exp.match(x^5)
    if str(m)!='None':
        return 'pow('+str(sage2C(m[x]))+',5)'
    m=exp.match(x^6)
    if str(m)!='None':
        return 'pow('+str(sage2C(m[x]))+',6)'
    m=exp.match(x^(3/2))
    if str(m)!='None':
        return 'pow(sqrt('+str(sage2C(m[x]))+'),3)'
    m=exp.match(x^(-3/2))
    if str(m)!='None':
        return '1./pow(sqrt('+str(sage2C(m[x]))+'),3)'
    m=exp.match(x^(-2))
    if str(m)!='None':
        return 'I(1)/pow('+str(sage2C(m[x]))+',2)'
    m=exp.match(arctan(x))    
    if str(m)!='None':
        return '  atan('+str(sage2C(m[x]))+')  ' 
    m=exp.match(sqrt(x))
    if str(m)!='None':
        return 'sqrt('+str(sage2C(m[x]))+')'
    print(exp)
    print ('!!!!!!!!!!!!!!!! ouuuuups!!!!')

# <> -> !=
var('ra','rb','rc','rd','ab','ac','ad','bc','bd','cd')
d = 8*(ra^3*
       arctan(
           sqrt((ac^2*bd^2 + (bc^2 - bd^2)*ad^2 - (ab^2 - ac^2 - ad^2 - bc^2 - bd^2 + cd^2)*cd^2)*ab^2 - (ab^2*bc^2 - (bc^2 - cd^2)*ad^2 + (ac^2 - ad^2 - bc^2 + bd^2 - cd^2)*bd^2)*ac^2 - (bd^2*cd^2 + (ad^2 + bc^2 - bd^2 - cd^2)*ad^2)*bc^2)/(2*ab*ac*ad + (ac^2 + ad^2 - cd^2)*ab + (ab^2 + ad^2 - bd^2)*ac + (ab^2 + ac^2 - bc^2)*ad)
       ) + rb^3*
       arctan(
           sqrt((ac^2*bd^2 - (ac^2 - cd^2)*ab^2 + (ab^2 + ac^2 - ad^2 - bc^2 + bd^2 + cd^2)*ad^2)*bc^2 - (bc^2*cd^2 + (ad^2 - cd^2)*ab^2 - (ab^2 - ac^2 + ad^2 - bd^2 + cd^2)*ac^2)*bd^2 - (ac^2*ad^2 + (ab^2 - ac^2 - ad^2 + cd^2)*ab^2)*cd^2)/(2*ab*bc*bd + (bc^2 + bd^2 - cd^2)*ab + (ab^2 - ad^2 + bd^2)*bc + (ab^2 - ac^2 + bc^2)*bd)
       )  + rc^3*
       arctan(
           sqrt(-(ad^2*cd^2 + (ab^2 - ad^2)*bc^2 - (ab^2 - ac^2 + ad^2 + bc^2 - bd^2)*bd^2)*ac^2 - (ab^2*bd^2 - (ab^2 - ad^2 - bc^2 + bd^2)*bc^2)*ad^2 + (ac^2*bd^2 - (ab^2 - ac^2 - ad^2 - bc^2 - bd^2 + cd^2)*ab^2 + (ad^2 - bd^2)*bc^2)*cd^2)/(2*ac*bc*cd + (bc^2 - bd^2 + cd^2)*ac + (ac^2 - ad^2 + cd^2)*bc - (ab^2 - ac^2 - bc^2)*cd))
       + rd^3*
       arctan(
           sqrt(-(ac^2*bc^2 + (ab^2 - ac^2 - bc^2 + cd^2)*cd^2)*ab^2 + (ac^2*bd^2 + (ab^2 + ac^2 - ad^2 - bc^2 + bd^2 + cd^2)*bc^2 + (ab^2 - ac^2)*cd^2)*ad^2 - (ab^2*ad^2 - (ab^2 - ac^2 + bc^2 - bd^2 + cd^2)*ac^2 - (ab^2 - bc^2)*cd^2)*bd^2)/(2*ad*bd*cd - (bc^2 - bd^2 - cd^2)*ad - (ac^2 - ad^2 - cd^2)*bd - (ab^2 - ad^2 - bd^2)*cd)
       ))/sqrt((ac^2*bd^2 + (bc^2 - bd^2)*ad^2 - (ab^2 - ac^2 - ad^2 - bc^2 - bd^2 + cd^2)*cd^2)*ab^2 - (ab^2*bc^2 - (bc^2 - cd^2)*ad^2 + (ac^2 - ad^2 - bc^2 + bd^2 - cd^2)*bd^2)*ac^2 - (bd^2*cd^2 + (ad^2 + bc^2 - bd^2 - cd^2)*ad^2)*bc^2)

den = (ac^2*bd^2 + (bc^2 - bd^2)*ad^2 - (ab^2 - ac^2 - ad^2 - bc^2 - bd^2 + cd^2)*cd^2)*ab^2 - (ab^2*bc^2 - (bc^2 - cd^2)*ad^2 + (ac^2 - ad^2 - bc^2 + bd^2 - cd^2)*bd^2)*ac^2 - (bd^2*cd^2 + (ad^2 + bc^2 - bd^2 - cd^2)*ad^2)*bc^2

atA = 'sqrt(' + sage2C((ac^2*bd^2 + (bc^2 - bd^2)*ad^2 - (ab^2 - ac^2 - ad^2 - bc^2 - bd^2 + cd^2)*cd^2)*ab^2 - (ab^2*bc^2 - (bc^2 - cd^2)*ad^2 + (ac^2 - ad^2 - bc^2 + bd^2 - cd^2)*bd^2)*ac^2 - (bd^2*cd^2 + (ad^2 + bc^2 - bd^2 - cd^2)*ad^2)*bc^2) +')' +'/(' + sage2C(2*ab*ac*ad + (ac^2 + ad^2 - cd^2)*ab + (ab^2 + ad^2 - bd^2)*ac + (ab^2 + ac^2 - bc^2)*ad) + ')'

atB = 'sqrt(' + sage2C((ac^2*bd^2 - (ac^2 - cd^2)*ab^2 + (ab^2 + ac^2 - ad^2 - bc^2 + bd^2 + cd^2)*ad^2)*bc^2 - (bc^2*cd^2 + (ad^2 - cd^2)*ab^2 - (ab^2 - ac^2 + ad^2 - bd^2 + cd^2)*ac^2)*bd^2 - (ac^2*ad^2 + (ab^2 - ac^2 - ad^2 + cd^2)*ab^2)*cd^2) +')' +'/(' + sage2C(2*ab*bc*bd + (bc^2 + bd^2 - cd^2)*ab + (ab^2 - ad^2 + bd^2)*bc + (ab^2 - ac^2 + bc^2)*bd) + ')'

atC = 'sqrt(' + sage2C(-(ad^2*cd^2 + (ab^2 - ad^2)*bc^2 - (ab^2 - ac^2 + ad^2 + bc^2 - bd^2)*bd^2)*ac^2 - (ab^2*bd^2 - (ab^2 - ad^2 - bc^2 + bd^2)*bc^2)*ad^2 + (ac^2*bd^2 - (ab^2 - ac^2 - ad^2 - bc^2 -
                                                                                                                                                                              bd^2 + cd^2)*ab^2 + (ad^2 - bd^2)*bc^2)*cd^2) +')' +'/(' + sage2C(2*ac*bc*cd + (bc^2 - bd^2 + cd^2)*ac + (ac^2 - ad^2 + cd^2)*bc - (ab^2 - ac^2 - bc^2)*cd) + ')'

atD = 'sqrt(' + sage2C(-(ac^2*bc^2 + (ab^2 - ac^2 - bc^2 + cd^2)*cd^2)*ab^2 + (ac^2*bd^2 + (ab^2 + ac^2 - ad^2 - bc^2 + bd^2 + cd^2)*bc^2 + (ab^2 - ac^2)*cd^2)*ad^2 - (ab^2*ad^2 - (ab^2 - ac^2 + bc^2 - bd^2 + cd^2)*ac^2 - (ab^2 - bc^2)*cd^2)*bd^2) +')' +'/(' + sage2C(2*ad*bd*cd - (bc^2 - bd^2 - cd^2)*ad - (ac^2 - ad^2 - cd^2)*bd - (ab^2 - ad^2 - bd^2)*cd) + ')'

D = '8*(' + 'pow(ra,3)*atan(' + atA + ')' + '+pow(rb,3)*atan(' + atB + ')' + '+pow(rc,3)*atan(' + atC + ')' + '+pow(rd,3)*atan(' + atD + ')' +')/'+ 'sqrt(' + sage2C(den) + ')'

square = lambda x: x**2
atam = lambda x: arctan(x)
I = lambda x: RIF(x)

#print(D)

var('r')

