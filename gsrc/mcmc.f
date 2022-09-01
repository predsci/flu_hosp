      subroutine mcmc(vecN, vecI0, nparam, nb, par,
     $     parmin, parmax, step, ilogflag, imask, iseed,
     $     wght, obs, gamaObs, trickle, vecTcalc, nTimes, ndays,
     $     rtn_daily, rec_daily, hos_daily, nMCMC, ithin, tab)

      
      implicit none

      integer nparam, nb, nparamtot
      integer nTimes, ndays
      real*8 par((nparam+2*nb))
      real*8 parmin((nparam+2*nb)), parmax((nparam+2*nb))
      real*8 step((nparam+2*nb))
      integer ilogflag((nparam+2*nb)), imask((nparam+2*nb)), iseed
      real*8 wght(ndays), obs(ndays), gamaObs(ndays)
      real*8 vecN, vecI0
      real*8 trickle
      real*8 vecTcalc(nTimes)

      real*8 sum_ts
      
!     Declare items to be returned
      
      real*8 rtn_daily(ndays)
      real*8 rec_daily(ndays)
      real*8 hos_daily(ndays)
      
!     mcmc items
      integer nMCMC, ithin
      real  tab(nMCMC/ithin,(nparam+2*nb+1))
      
      real*8 obsLLK, curLLK, curMin, newLLK
      real*8 calcfit
      
      integer ii
      integer nopt, iopt((nparam+2*nb))
      real*8 savepar((nparam+2*nb)), curpars((nparam+2*nb))
      real*8 copypar((nparam+2*nb)), parupdt((nparam+2*nb))
      real*8 savestep((nparam+2*nb))
      real*8 hosBest(ndays), parBest((nparam+2*nb))
      
      integer nlines, ionep, iadapt, iaccept, icount
      real*8 scale, range_min, range_max, step_min, step_max, myaccept
      real*8 rv
      logical accept

!     Total number of parameters

      nparamtot = nparam + nb * 2
      
! Number of lines in the tab array

      nlines = nMCMC/ithin

      tab = 0.0
! For an adaptive size MCMC - decide if we need to update step size every 1% 
! or 1,000 steps
      scale = 2.0d0
      ionep = int(nMCMC * 0.01)
      ionep = min(1000, ionep)
      ionep = max(ionep, 1)
      range_min = 0.20d0
      range_max = 0.30d0
      step_max = 1.d0
      step_min = 1e-5
      iadapt = 0

      
! ran1 works with -iseed

      iseed = -abs(iseed)

      call srand(iseed)
      
!     book keeping - who gets optimized and who does not

      nopt = 0

      do ii = 1, nparamtot
         if (imask(ii) > 0) then
            nopt = nopt + 1
            iopt(nopt) = ii
         endif
      enddo
      

!     ensure that the last time change does not exceed number of days
      par(nparamtot) = 0.0d0
      sum_ts = sum(par((nparam+nb+1):nparamtot))

      if (sum_ts > ndays) Then
         par((nparam+nb+1):nparamtot) =
     $        par((nparam+nb+1):nparamtot) * dble(ndays)/sum_ts
      endif      
!     calculate initial profile      
      curpars = par
      
      call detsirh(vecN, vecI0, nparam, nb, par,
     $     trickle, vecTcalc, nTimes, ndays,
     $     rtn_daily, rec_daily,hos_daily)

!
! calculate Likelihood of the solution
!

!     This is the best scenario basically
      obsLLK = calcfit(obs,gamaObs,obs+1e-10,wght,ndays)
!     Initial LLK     
      curLLK = calcfit(obs,gamaObs,hos_daily,wght,ndays)

!
! MCMC loop starts here 
!
      curMIN = curLLK
      hosBest = hos_daily 

      iaccept = 0
      icount = 0

      savestep = step
     
      parupdt = curpars
      
      savepar = curpars
      copypar = curpars

      do ii = 1, nMCMC
         
         savepar = curpars
         copypar = curpars       
!     Propose a new value for each parameter we are optimizing 
!     half the steps will be small and half will be the given size
         call fnProposeParamUpdates(nparamtot,copypar,
     $        parmin,parmax,step,
     $        ilogflag,parupdt,iseed,nopt,iopt(1:nopt))

         parupdt(nparamtot) = 0.0d0
         sum_ts = sum(parupdt((nparam+nb+1):nparamtot))

         if (sum_ts > ndays) Then
            parupdt((nparam+nb+1):nparamtot) =
     $           parupdt((nparam+nb+1):nparamtot) * dble(ndays)/sum_ts
         endif
         
         curpars(iopt(1:nopt)) = parupdt(iopt(1:nopt))

!     calculate a new profile 
           call detsirh(vecN, vecI0, nparam, nb, curpars,
     $        trickle, vecTcalc, nTimes, ndays,
     $        rtn_daily, rec_daily,hos_daily)
         
           newLLK = calcfit(obs,gamaObs,hos_daily,wght,ndays)

           rv = rand()

           if (exp((curLLK - newLLK)) .gt. rv ) Then
              accept = .true.
              iaccept = iaccept + 1
              iadapt = iadapt + 1
              curLLK = newLLK
              savepar = curpars
              hosBest = hos_daily
              parbest = curpars
           else
              curpars = savepar
           endif

        if (mod(ii,ithin) .eq. 0) Then
            icount = icount + 1
            tab(icount,1:nparamtot)   = real(curpars)
            tab(icount,(nparamtot+1)) = real(curLLK)
           
         endif

         if (mod(ii, 1000) .eq. 0) then
            print*, ii, real(curLLK),
     $            real(curpars((nparam+1):(nparam+2*nb-1)))
         end if
         
         if (mod(ii,ionep) .eq. 0) Then
            
            myaccept = dble(iadapt)/dble(ionep)
            if (myaccept > range_max .and.
     $           all(step(iopt(1:nopt))*scale < step_max)) Then
                 step(iopt(1:nopt))= step(iopt(1:nopt))*scale
               endif
               if (myaccept < range_min .and.
     $         all(step(iopt(1:nopt))/scale > step_min)) Then
                 step(iopt(1:nopt))= step(iopt(1:nopt))/scale
               endif
               iadapt = 0    !reset acceptance number

         endif
         

      enddo

      par = parBest  ! return the best estimate for the parameters
      return
      end subroutine mcmc

!--------------------------------------------------------------------------------
        function calcfit(y,gamay,x,wght,ndata)

        implicit none

        integer ndata, i
        real*8 y(ndata),x(ndata), gamay(ndata)
        real*8 wght(ndata)
        real*8 xi, yi,sum,val,calcfit


c x is the simulated data
c y is the base profile

C calculate the P(yi,xi)


        sum = 0.0
        do i=1,ndata
           
           yi = y(i)
           xi = x(i)
           
           val = yi * log(xi) - xi - gamay(i)
           sum = sum + val  * wght(i) 
           
        enddo
        sum = -sum   
        
        calcfit = sum 
           
        return
        end function calcfit


!----------------------------------------------------------------

      subroutine fnProposeParamUpdates(nparam,curval,
     $     valmin,valmax,step,ilogflag,
     $     parupdt,iseed,nopt,iopt)

      implicit none
      integer nparam, nopt, iopt(nopt)
      real*8 curval(nparam),valmin(nparam),valmax(nparam)
      real*8 parupdt(nparam),step(nparam)
      real*8 x, rv, rtn
      real*8 ran1
      real*8 SR_to_unit, SR_from_unit
      integer iseed
      integer i,j
      
      integer ilogflag(nparam)
      external SR_to_unit, SR_from_unit

       do j = 1, nopt
          i = iopt(j)
          parupdt(i)=curval(i)

          rv = rand()

          rv = (rv - 0.50d0)*step(i)

! convert to a zero - one scale

          x = SR_to_unit(curval(i),valmin(i),valmax(i),
     $         ilogflag(i))

          x = x + rv

      if (x .lt. 0.0d0) x = 1.0d0 + x
      if (x .gt. 1.0d0) x = x - 1.0d0

! Do not use period boundary conditions here but rather re-smaple
c$$$      if (x .le. 0.0d0 .or. x .ge. 1.0d0) go to 101


! bring value back to original scale
      
         rtn = SR_from_unit(x,valmin(i),valmax(i),
     $        ilogflag(i))

         parupdt(i) = rtn

      enddo

      return
      end subroutine fnProposeParamUpdates


c----------------------------------------------------------------

      function SR_to_unit(y,ymin,ymax, ilogflag)

      implicit none
      real*8 y, ymin,ymax,rtn,SR_to_Unit
      integer ilogflag

      if (ilogflag == 1) Then
         rtn = (log10(y) - log10(ymin)) /
     $        (log10(ymax)-log10(ymin))
      else
         rtn = (y - ymin)/(ymax - ymin)
      endif

      SR_to_unit = rtn
 
      return
      end function SR_to_unit

c----------------------------------------------------------------

      function SR_from_unit(x,ymin,ymax, ilogflag)
      
      implicit none
   
      real*8 x,ymin,ymax,rtn,SR_from_unit
      integer ilogflag


      if (ilogflag == 1) Then
         rtn = ymin * 
     $        10.0**(x*(log10(ymax)-log10(ymin)))

      else
         rtn = ymin + (ymax-ymin)*x
      endif

      SR_from_unit = rtn

      return
      end function SR_from_unit

C--------------------------------------------------------------------------------
      
