#!/usr/bin/env python
#-*- coding:utf-8 -*-

import math

# Global variable equivalent to swift.inc
NPLMAX = 51   # max number of planets, including the Sun
NTPMAX = 1001 # max number of test particles
    
#Size of the test particle integer status flag
#NSTATP Number of status parameters
NSTATP = 3
#NSTAT Number of status parameters
NSTAT = NSTATP + NPLMAX - 1 # include one for @ planet  

# Size of the test particle integer status flag
# integer NSTATR    
NSTATR = NSTAT # io_init_tp assumes NSTAT==NSTATR

# Convergence criteria for danby
# real DANBYAC , DANBYB
DANBYAC = 1.0e-14
DANBYB = 1.0e-13

# loop limits in the Laguerre attempts
# integer NLAG1, NLAG2
NLAG1 = 50
NLAG2 = 400

# A small number
# real*8 TINY
TINY = 4.0e-15

# trig stuff
# real*8 PI,TWOPI,PIBY2,DEGRAD
PI = 3.14159265358979e0
TWOPI = 2.0e0 * PI
PIBY2 = PI/2.0e0
DEGRAD = 180.0e0 / PI

def orbel_el2xv(gm, ialpha, a, e, inc, capom, omega, capm):
    """
    *****************************************************************************
    *                          ORBEL_EL2XV.F
    *****************************************************************************
    *     PURPOSE: To compute cartesian positions and velocities given
    *               central mass, ialpha ( = +1 for hyp., 0 for para. and
    *               -1 for ellipse), and orbital elements.
    C       input:
    c            gm       ==> G times central mass (real scalar)
    c            ialpha   ==> conic section type ( see PURPOSE, integer scalar)
    C            a        ==> semi-major axis or pericentric distance if a parabola
    c                          (real scalar)
    c            e        ==> eccentricity (real scalar)
    C            inc      ==> inclination  (real scalar)
    C            capom    ==> longitude of ascending node (real scalar)
    C            omega    ==> argument of perihelion (real scalar)
    C            capm     ==> mean anomoly(real scalar)
    *       
    c       Output:
    c            x,y,z    ==>  position of object (real scalars)
    c            vx,vy,vz ==>  velocity of object (real scalars)
    c
    *     ALGORITHM:  See Fitzpatrick "Principles of Cel. Mech."
    *     REMARKS: All angles are in RADIANS
    *       
    *     AUTHOR:  M. Duncan.
    *     DATE WRITTEN:  May 11, 1992.
    *     REVISIONS: May 26 - now use better Kepler solver for ellipses
    *                 and hyperbolae called EHYBRID.F and FHYBRID.F
    ***********************************************************************
    """

    if e < 0.0:
        print(" ERROR in orbel_el2xv: e<0, setting e=0!")
        e = 0.0

    # Check for inconsistencies between ialpha and e
    em1 = e - 1.0e0
    if ((ialpha == 0) and (abs(em1) > TINY))  or \
        ((ialpha < 0) and (e > 1.0e0)) or \
        ((ialpha > 0) and (e < 1.0e0)):
        print("ERROR in orbel_el2xv: ialpha and e inconsistent")
        print("ialpha = ", ialpha)
        print("e = ", e)

    # Generate rotation matrices (on p. 42 of Fitzpatrick)
    sp, cp = orbel_scget(omega)
    so, co = orbel_scget(capom)
    si, ci = orbel_scget(inc)
    
    d11 = cp*co - sp*so*ci
    d12 = cp*so + sp*co*ci
    d13 = sp*si
    d21 = -sp*co - cp*so*ci
    d22 = -sp*so + cp*co*ci
    d23 = cp*si

    # Get the other quantities depending on orbit type ( i.e. IALPHA)
    if ialpha == -1:
        cape = orbel_ehybrid(e,capm)
        scap, ccap = orbel_scget(cape)
        sqe = math.sqrt(1.0e0 - e*e)
        sqgma = math.sqrt(gm * a)
        xfac1 = a*(ccap - e)
        xfac2 = a*sqe*scap
        ri = 1.0e0 / (a*(1.0e0 - e*ccap))
        vfac1 = -ri * sqgma * scap
        vfac2 = ri * sqgma * sqe * ccap

    elif ialpha == 1:
        capf = orbel_fhybrid(e,capm)
        shcap, chcap = orbel_schget(capf)
        sqe = math.sqrt(e*e - 1.0e0 )
        sqgma = math.sqrt(gm*a)
        xfac1 = a*(e - chcap)
        xfac2 = a*sqe*shcap
        ri = 1.0e0 / (a*(e*chcap - 1.0e0))
        vfac1 = -ri * sqgma * shcap
        vfac2 = ri * sqgma * sqe * chcap

    else:
        zpara = orbel_zget(capm)
        sqgma = math.sqrt(2.0e0*gm*a)
        xfac1 = a*(1.0e0 - zpara*zpara)
        xfac2 = 2.0e0*a*zpara
        ri = 1.0e0/(a*(1.0e0 + zpara*zpara))
        vfac1 = -ri * sqgma * zpara
        vfac2 = ri * sqgma

    x =  d11*xfac1 + d21*xfac2
    y =  d12*xfac1 + d22*xfac2
    z =  d13*xfac1 + d23*xfac2
    vx = d11*vfac1 + d21*vfac2
    vy = d12*vfac1 + d22*vfac2
    vz = d13*vfac1 + d23*vfac2

    return x, y, z, vx, vy, vz

    # End orbel_el2xv

def orbel_scget(angle):
    """
    ***********************************************************************
    c                         ORBEL_SCGET.F
    ***********************************************************************
    *     PURPOSE:  Given an angle, efficiently compute sin and cos.
    *
    *        Input:
    *             angle ==> angle in radians (real scalar)
    *        
    *        Output:
    *             sx    ==>  sin(angle)  (real scalar)
    *             cx    ==>  cos(angle)  (real scalar)
    *
    *     ALGORITHM: Obvious from the code 
    *     REMARKS: The HP 700 series won't return correct answers for sin
    *       and cos if the angle is bigger than 3e7. We first reduce it
    *       to the range [0,2pi) and use the sqrt rather than cos (it's faster)
    *       BE SURE THE ANGLE IS IN RADIANS - NOT DEGREES!
    *     AUTHOR:  M. Duncan.
    *     DATE WRITTEN:  May 6, 1992.
    *     REVISIONS: 
    ***********************************************************************
    """

    PI3BY2 = 1.5e0*PI
    nper = int(angle / TWOPI)
    x = angle - nper * TWOPI
    if x < 0.0e0:
        x = x + TWOPI
    sx = math.sin(x)
    cx = math.sqrt(1.0e0 - sx*sx)
    if x > PIBY2 and x < PI3BY2:
        cx = -cx
    
    return sx, cx
    # end orbel_scget

def orbel_ehybrid(e,m):
    """
    ***********************************************************************
    c                    ORBEL_EHYBRID.F
    ***********************************************************************
    *     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
    *
    *             Input:
    *                           e ==> eccentricity anomaly. (real scalar)
    *                           m ==> mean anomaly. (real scalar)
    *             Returns:
    *              orbel_ehybrid ==>  eccentric anomaly. (real scalar)
    *
    *     ALGORITHM: For e < 0.18 uses fast routine ESOLMD 
    *                For larger e but less than 0.8, uses EGET
    *                For e > 0.8 uses EHIE
    *     REMARKS: Only EHIE brings M and E into range (0,TWOPI)
    *     AUTHOR: M. Duncan 
    *     DATE WRITTEN: May 25,1992.
    *     REVISIONS: 2/26/93 hfl
    ***********************************************************************
    """
    if e < 0.18e0:
        orbel_ehybrid = orbel_esolmd(e,m)
    else:
        if e <= 0.8e0:
            orbel_ehybrid = orbel_eget(e,m)
        else:
            orbel_ehybrid = orbel_ehie(e,m) 

    return orbel_ehybrid

    # End orbel_ehybrid

def orbel_fhybrid(e, n):
    """
    ***********************************************************************
    c                    ORBEL_FHYBRID.F
    ***********************************************************************
    *     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
    *
    *             Input:
    *                           e ==> eccentricity anomaly. (real scalar)
    *                           n ==> hyperbola mean anomaly. (real scalar)
    *             Returns:
    *               orbel_fhybrid ==>  eccentric anomaly. (real scalar)
    *
    *     ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON 
    *                For larger N, uses FGET
    *     REMARKS: 
    *     AUTHOR: M. Duncan 
    *     DATE WRITTEN: May 26,1992.
    *     REVISIONS: 
    *     REVISIONS: 2/26/93 hfl
    ***********************************************************************
    """
    
    abn = n
    
    if n < 0.0e0:
        abn = -abn

    if abn < 0.636e0*e -0.6e0:
        orbel_fhybrid = orbel_flon(e,n)
    else:
        orbel_fhybrid = orbel_fget(e,n)

    return orbel_fhybrid

    # end orbel_fhybrid

def orbel_zget(q):
    """
    ***********************************************************************
    c                    ORBEL_ZGET.F
    ***********************************************************************
    *     PURPOSE:  Solves the equivalent of Kepler's eqn. for a parabola 
    *          given Q (Fitz. notation.)
    *
    *             Input:
    *                           q ==>  parabola mean anomaly. (real scalar)
    *             Returns:
    *                  orbel_zget ==>  eccentric anomaly. (real scalar)
    *
    *     ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
    *     REMARKS: For a parabola we can solve analytically.
    *     AUTHOR: M. Duncan 
    *     DATE WRITTEN: May 11, 1992.
    *     REVISIONS: May 27 - corrected it for negative Q and use power
    *             series for small Q.
    ***********************************************************************
    """
    iflag = 0
    
    if q < 0.0e0:
        iflag = 1
        q = -q
    
    if q < 1.0e-3:
        orbel_zget = q*(1.0e0 - (q*q/3.0e0)*(1.0e0 -q*q))
    else:
        x = 0.5e0*(3.0e0*q + sqrt(9.0e0*(q**2) +4.0e0))
        tmp = x**(1.0e0 / 3.0e0)
        orbel_zget = tmp - 1.0e0 / tmp
    
    if iflag == 1:
        orbel_zget = -orbel_zget
        q = -q
       
    return orbel_zget
    # end orbel_zget

def orbel_esolmd(e,m):
    """
    ***********************************************************************
    c                    ORBEL_ESOLMD.F
    ***********************************************************************
    *     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
    *
    *             Input:
    *                           e ==> eccentricity anomaly. (real scalar)
    *                           m ==> mean anomaly. (real scalar)
    *             Returns:
    *                orbel_esolmd ==>  eccentric anomaly. (real scalar)
    *
    *     ALGORITHM: Some sort of quartic convergence from Wisdom. 
    *     REMARKS: ONLY GOOD FOR SMALL ECCENTRICITY SINCE IT ONLY
    *         ITERATES ONCE. (GOOD FOR PLANET CALCS.)
    *               ALSO DOES NOT PUT M OR E BETWEEN 0. AND 2*PI 
    *     INCLUDES: needs SCGET.F
    *     AUTHOR: M. Duncan 
    *     DATE WRITTEN: May 7, 1992.
    *     REVISIONS: 2/26/93 hfl
    ***********************************************************************
    """
    sm, cm = orbel_scget(m)
    x = m + e*sm*( 1.0e0 + e*( cm + e*( 1.0e0 -1.5e0*sm*sm)))

    sx, cx = orbel_scget(x)
    es = e*sx
    ec = e*cx
    f = x - es  - m
    fp = 1.0e0 - ec 
    fpp = es 
    fppp = ec 
    dx = -f/fp
    dx = -f/(fp + dx*fpp / 2.0e0)
    dx = -f/(fp + dx*fpp / 2.0e0 + dx*dx*fppp / 6.0e0)

    orbel_esolmd = x + dx

    return orbel_esolmd
    # end orbel_esolmd

def orbel_eget(e,m):
    """
    c--------------------------------------------------------------------
    ***********************************************************************
    c                    ORBEL_EGET.F
    ***********************************************************************
    *     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
    *
    *             Input:
    *                           e ==> eccentricity anomaly. (real scalar)
    *                           m ==> mean anomaly. (real scalar)
    *             Returns:
    *                  orbel_eget ==>  eccentric anomaly. (real scalar)
    *
    *     ALGORITHM: Quartic convergence from Danby
    *     REMARKS: For results very near roundoff, give it M between
    *           0 and 2*pi. One can condition M before calling EGET
    *           by calling my double precision function MOD2PI(M). 
    *           This is not done within the routine to speed it up
    *           and because it works fine even for large M.
    *     AUTHOR: M. Duncan 
    *     DATE WRITTEN: May 7, 1992.
    *     REVISIONS: May 21, 1992.  Now have it go through EXACTLY two iterations
    *                with the premise that it will only be called if
    *                we have an ellipse with e between 0.15 and 0.8
    ***********************************************************************
    """

    sm, cm = orbel_scget(m)

    # begin with a guess accurate to order ecc**3       
    x = m + e*sm*( 1.0e0 + e*( cm + e*( 1.0e0 -1.5e0*sm*sm)))

    # Go through one iteration for improved estimate
    sx, cx = orbel_scget(x)
    es = e*sx
    ec = e*cx
    f = x - es  - m
    fp = 1.0e0 - ec 
    fpp = es 
    fppp = ec 
    dx = -f/fp
    dx = -f/(fp + dx*fpp/2.0e0)
    dx = -f/(fp + dx*fpp/2.0e0 + dx*dx*fppp/6.0e0)
    orbel_eget = x + dx

    # c Do another iteration.
    # c For m between 0 and 2*pi this seems to be enough to
    # c get near roundoff error for eccentricities between 0 and 0.8

    x = orbel_eget
    sx, cx = orbel_scget(x)
    es = e*sx
    ec = e*cx
    f = x - es  - m
    fp = 1.0e0 - ec 
    fpp = es 
    fppp = ec 
    dx = -f / fp
    dx = -f / (fp + dx*fpp / 2.0e0)
    dx = -f / (fp + dx*fpp / 2.0e0 + dx*dx*fppp / 6.0e0)

    orbel_eget = x + dx

    return orbel_eget
    # end  ! orbel_eget

def orbel_ehie(e,m):
    """
    c---------------------------------------------------------------------
    ***********************************************************************
    c                    ORBEL_EHIE.F
    ***********************************************************************
    *     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
    *
    *             Input:
    *                           e ==> eccentricity anomaly. (real scalar)
    *                           m ==> mean anomaly. (real scalar)
    *             Returns:
    *              orbel_ehybrid ==>  eccentric anomaly. (real scalar)
    *
    *     ALGORITHM: Use Danby's quartic for 3 iterations. 
    *                Eqn. is f(x) = x - e*sin(x+M). Note  that
    *                E = x + M. First guess is very good for e near 1.
    *                Need to first get M between 0. and PI and use
    *               symmetry to return right answer if M between PI and 2PI
    *     REMARKS: Modifies M so that both E and M are in range (0,TWOPI)
    *     AUTHOR: M. Duncan 
    *     DATE WRITTEN: May 25,1992.
    *     REVISIONS: 
    ***********************************************************************
    """

    NMAX = 3

    # In this section, bring M into the range (0,TWOPI) and if
    # the result is greater than PI, solve for (TWOPI - M).
    iflag = 0
    nper = int(m / TWOPI)
    m = m - nper*TWOPI
    if m < 0.0e0:
        m = m + TWOPI

    if m > PI:
        m = TWOPI - m
        iflag = 1

    # Make a first guess that works well for e near 1.
    x = (6.0e0*m)**(1.0e0/3.0e0) - m
    niter = 0

    # Iteration loop
    for niter in range(1, NMAX):
        sa, ca = orbel_scget(x + m)
        esa = e*sa
        eca = e*ca
        f = x - esa
        fp = 1.0e0 -eca
        dx = -f/fp
        dx = -f/(fp + 0.50e0*dx*esa)
        dx = -f/(fp + 0.50e0*dx*(esa+0.3333333333333333e0*eca*dx))
        x = x + dx

    orbel_ehie = m + x

    if iflag == 1:
        orbel_ehie = TWOPI - orbel_ehie
        m = TWOPI - m

    return orbel_ehie
    # end orbel_ehie

def orbel_flon(e,capn):
    """
    ***********************************************************************
    c                    ORBEL_FLON.F
    ***********************************************************************
    *     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
    *
    *             Input:
    *                           e ==> eccentricity anomaly. (real scalar)
    *                        capn ==> hyperbola mean anomaly. (real scalar)
    *             Returns:
    *                  orbel_flon ==>  eccentric anomaly. (real scalar)
    *
    *     ALGORITHM: Uses power series for N in terms of F and Newton,s method
    *     REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)
    *     AUTHOR: M. Duncan 
    *     DATE WRITTEN: May 26, 1992.
    *     REVISIONS: 
    ***********************************************************************
    """
    IMAX = 10
    a11 = 156.0e0
    a9 = 17160.0e0
    a7 = 1235520.0e0
    a5 = 51891840.0e0
    a3 = 1037836800.0e0
    b11 = 11.0e0*a11
    b9 = 9.0e0*a9
    b7 = 7.0e0*a7
    b5 = 5.0e0*a5
    b3 = 3.0e0*a3


    # Function to solve "Kepler's eqn" for F (here called
    # x) for given e and CAPN. Only good for smallish CAPN 

    iflag = 0
    if capn < 0.0e0:
        iflag = 1
        capn = -capn

    a1 = 6227020800.0e0 * (1.0e0 - 1.0e0/e)
    a0 = -6227020800.0e0*capn / e
    b1 = a1

    #  Set iflag nonzero if capn < 0., in which case solve for -capn
    #  and change the sign of the final answer for F.
    #  Begin with a reasonable guess based on solving the cubic for small F       


    a = 6.0e0*(e-1.0e0)/e
    b = -6.0e0*capn/e
    sq = math.sqrt(0.25*b*b +a*a*a/27.0e0)
    biga = (-0.5*b + sq)**0.3333333333333333e0
    bigb = -(+0.5*b + sq)**0.3333333333333333e0
    x = biga + bigb
    # write(6,*) 'cubic = ',x**3 +a*x +b
    orbel_flon = x
    # If capn is tiny (or zero) no need to go further than cubic even for
    # e =1.
    if capn < TINY:
        if iflag == 1:
            orbel_flon = -orbel_flon
            capn = -capn
            return orbel_flon
    for i in range(1, IMAX):
        x2 = x*x
        f = a0 +x*(a1+x2*(a3+x2*(a5+x2*(a7+x2*(a9+x2*(a11+x2))))))
        fp = b1 +x2*(b3+x2*(b5+x2*(b7+x2*(b9+x2*(b11 + 13.0e0*x2)))))   
        dx = -f/fp
    # write(6,*) 'i,dx,x,f : '
    # write(6,432) i,dx,x,f
    # 432         format(1x,i3,3(2x,1p1e22.15))
    orbel_flon = x + dx
    # If we have converged here there's no point in going on
    if abs(dx) <= TINY:
        if iflag == 1:
            orbel_flon = -orbel_flon
            capn = -capn
            return orbel_flon
    x = orbel_flon     

    #  Abnormal return here - we've gone thru the loop 
    # IMAX times without convergence
    if iflag == 1:
        orbel_flon = -orbel_flon
        capn = -capn
    print("FLON: RETURNING WITHOUT COMPLETE CONVERGENCE")   
    diff = e*sinh(orbel_flon) - orbel_flon - capn
    print("N, F, ecc*sinh(F) - F - N : ")
    print(capn,orbel_flon,diff)
    return orbel_flon
    # end orbel_flon

def orbel_fget(e,capn):
    """
    ***********************************************************************
    c                    ORBEL_FGET.F
    ***********************************************************************
    *     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
    *
    *             Input:
    *                           e ==> eccentricity anomaly. (real scalar)
    *                        capn ==> hyperbola mean anomaly. (real scalar)
    *             Returns:
    *                  orbel_fget ==>  eccentric anomaly. (real scalar)
    *
    *     ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
    *           Cel. Mech. ".  Quartic convergence from Danby's book.
    *     REMARKS: 
    *     AUTHOR: M. Duncan 
    *     DATE WRITTEN: May 11, 1992.
    *     REVISIONS: 2/26/93 hfl
    ***********************************************************************
    """

    IMAX = 10
    if capn < 0.0e0:
        tmp = -2.0e0*capn/e + 1.80e0
        x = -log(tmp)
    else:
        tmp = +2.0e0*capn/e + 1.80e0
        x = math.log(tmp)
       
    orbel_fget = x

    for i in range(1, IMAX):
        shx, chx = orbel_schget(x)
        esh = e*shx
        ech = e*chx
        f = esh - x - capn
        # write(6,*) 'i,x,f : ',i,x,f
        fp = ech - 1.0e0  
        fpp = esh 
        fppp = ech 
        dx = -f/fp
        dx = -f/(fp + dx*fpp/2.0e0)
        dx = -f/(fp + dx*fpp/2.0e0 + dx*dx*fppp/6.0e0)
        orbel_fget = x + dx
        # If we have converged here there's no point in going on
        if abs(dx) <= TINY:
            return orbel_fget

        x = orbel_fget    

    print("FGET : RETURNING WITHOUT COMPLETE CONVERGENCE")
    return orbel_fget

    # end orbel_fget

def orbel_schget(angle):
    """
    ***********************************************************************
    c                         ORBEL_SCHGET.F
    ***********************************************************************
    *     PURPOSE:  Given an angle, efficiently compute sinh and cosh.
    *
    *        Input:
    *             angle ==> angle in radians (real scalar)
    *        
    *        Output:
    *             shx    ==>  sinh(angle)  (real scalar)
    *             chx    ==>  cosh(angle)  (real scalar)
    *
    *     ALGORITHM: Obvious from the code 
    *     REMARKS: Based on the routine SCGET for sine's and cosine's.
    *       We use the sqrt rather than cosh (it's faster)
    *       BE SURE THE ANGLE IS IN RADIANS AND IT CAN'T BE LARGER THAN 300
    *       OR OVERFLOWS WILL OCCUR!
    *     AUTHOR:  M. Duncan.
    *     DATE WRITTEN:  May 6, 1992.
    *     REVISIONS: 
    ***********************************************************************
    """
    shx = math.sinh(angle)
    chx = math.sqrt(1.0e0 + shx*shx)

    return shx, chx
    # end orbel_schget
