      subroutine pwr2stg(m,m1,n1,n2,x1,x2,muvec,nmu,muloc,rho0,rho1,
     1     nsel,nsim,pow,tmp,ord1,rnk1)
      integer m,m1,n1,n2,nmu,muloc(nmu),nsel,nsim,ord1(m),rnk1(m)
      double precision x1(m),x2(m),muvec(m),rho0,rho1,pow,tmp(m)

      double precision sqrho0,sqrho1,sqn1,sqn2,x0,yi,x2min,dnorm
      integer n,tsel,npow,ns,i

      call rndstart()
      n = n1 + n2
      sqn = sqrt(dfloat(n))
      sqn1 = sqrt(dfloat(n1))
      sqn2 = sqrt(dfloat(n2))
      sqrho0 = sqrt(1-rho0**2)
      sqrho1 = sqrt(1-rho1**2)

      npow = 0
      do 100 ns = 1,nsim
         x0 = dnorm()
         yi = dnorm()
         x1(1) = rho0*x0 + sqrho0*(yi + muvec(1)*sqn1)
         tmp(1) = x1(1)
         ord1(1) = 1
         do 10 i = 1,m-1
            yi = rho1*yi + sqrho1*dnorm()
            x1(i+1) = rho0*x0 + sqrho0*(yi + muvec(i+1)*sqn1)
            tmp(i+1) = x1(i+1)
            ord1(i+1) = i+1
 10      continue
         call qsort4(tmp,ord1,1,m)
         do 15 i = 1,m
            rnk1(ord1(i)) = m+1-i
 15      continue
         tsel = 0
         do 20 i = 1,nmu
            if (rnk1(muloc(i)).le.m1) tsel = tsel + 1
 20      continue
         if (tsel.lt.nsel) go to 100
         x0 = dnorm()
         yi = dnorm()
         x2(1) = (sqn1*x1(1) + sqn2*(rho0*x0 + sqrho0*(yi + 
     1        muvec(1)*sqn2)))/sqn
         x2min = x2(1)
         do 30 i = 1,m-1
            yi = rho1*yi + sqrho1*dnorm()
            x2(i+1) = (sqn1*x1(i+1) + sqn2*(rho0*x0 + sqrho0*(yi + 
     1           muvec(i+1)*sqn2)))/sqn
            x2min = min(x2min,x2(i+1))
 30      continue
         do 35 i = 1,m-m1
            x2(ord1(i)) = x2(ord1(i)) + x2min
 35      continue
         do 36 i=1,m
            tmp(i) = x2(i)
            ord1(i) = i
 36      continue
         call qsort4(tmp,ord1,1,m)
         do 40 i = 1,m
            rnk1(ord1(i)) = m+1-i
 40      continue
         tsel = 0
         do 45 i = 1,nmu
            if (rnk1(muloc(i)).le.nsel) tsel = tsel + 1
 45      continue
         if (tsel.eq.nsel) npow = npow + 1
 100  continue
      pow = dfloat(npow)/dfloat(nsim)
      call rndend()

      return
      end
