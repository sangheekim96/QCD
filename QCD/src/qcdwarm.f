      subroutine qcdupdate(v, x, n, tau, lambda, vout)

      integer i, j, n, stopflag
      double precision  v(1:n), x(1:n), vout, tau, lambda
      double precision  tmp, S1, Sx

      S1 = -lambda
      Sx = lambda
      vout = 0
      
      do 10 i=1,n
         if (x(i) .GT. 0) then
            S1 = S1 - x(i)*tau
         else
            S1 = S1 + x(i)*(1-tau)
         endif

         if (v(i) .LT. 0) then
            Sx = Sx+abs(x(i))
         endif
         
 10   continue

      if (abs(S1+Sx) .LT. lambda) goto 150
      
      stopflag = 0
      i = 1

 100  continue
      
         do 200 j=i+1,n
            if (v(i) .GT. v(j)) then
               tmp = v(i)
               v(i) = v(j)
               v(j) = tmp

               tmp = x(i)
               x(i) = x(j)
               x(j) = tmp
            endif
 200     continue

         tmp = S1 + abs(x(i))
         if (v(i) .EQ. 0) then
            tmp = tmp + 2*lambda
         endif
         
         if ((S1 .LT. 0) .AND. (tmp .GE. 0)) then
            vout = v(i)
            stopflag = 1
         endif
         S1 = tmp
         
         i = i+1
      if ((stopflag. EQ. 0) .AND. (i .LT. n)) goto 100

         
 150  return
      end


      subroutine qcdwarm(x, y, n, p, tau, lambda, b0, bhat, it, mt, tl)
      integer n, p
      integer i, j, k, it, mt
      double precision x(1:n, 1:p), xj(1:(n+1)), y(1:n), b0(1:p)
      double precision tau, lambda, bhat(1:p), r(1:n), rj(1:n)
      double precision vout, ro(1:(n+1)), tl, bdiff


      do 33 j=1,p
         bhat(j) = b0(j)
 33   continue

      it = 0
      

 77   continue
      
      it = it+1
      bdiff = 0.0
      
         do 88 i = 1,n
            r(i) = y(i)
            do 99 j=1,p
               r(i) = r(i)-x(i,j)*bhat(j)
 99         continue
 88      continue

         do 66 j=1,p
            
            do 55 i=1,n
               rj(i) = r(i)+x(i,j)*bhat(j)
               ro(i) = rj(i)/x(i,j)
               xj(i) = x(i,j)
 55         continue
            ro(n+1) = 0
            xj(n+1) = 0

            k = n+1
            vout = 0
            
            call qcdupdate(ro, xj, k, tau, lambda, vout)
            bdiff = bdiff+abs(bhat(j)-vout)
            bhat(j) = vout
            
            do 44 i=1,n
               r(i) = rj(i)-x(i,j)*bhat(j)
 44         continue

 66      continue
         
      if ((bdiff .GT. tl) .AND. (it .LE. mt)) goto 77
         
      
      return
      end
      
            
      

      
