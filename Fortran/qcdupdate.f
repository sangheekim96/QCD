      subroutine qcdupdate(v, x, n, tau, lambda, vout)

      integer i, j, n, stopflag
      double precision  v(1:n), x(1:n), vout, tau, lambda
      double precision  tmp, S1

      S1 = -lambda
      
      do 10 i=1,n
         if (x(i) .GT. 0) then
            S1 = S1 - x(i)*tau
         else
            S1 = S1 + x(i)*(1-tau)
         endif
 10   continue

      stopflag = 0
      i = 1
c     100  continue
      
      do 100 i=1,(n-1)
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
         
c         i = i+1
c      if ((stopflag. EQ. 0) .AND. (i .LT. n)) goto 100

 100     continue
         
      return
      end
      
            
      

      
