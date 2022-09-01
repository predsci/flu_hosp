      subroutine stosirh(vecN, vecI0, nparam, nb, par_tab,
     $     trickle, vecTcalc, nReals, nTimes, ndays,
     $     rtn_daily, rec_daily,hos_daily)

      implicit none
      
      integer nparam, nb
      integer nTimes, nReals, ndays

      real*8 wl
      parameter (wl = 3.0d0)
      
      real*8 par_tab(nReals, (nparam+nb*2))
      real*8 par(nparam+(nb*2))
      integer vecN, vecI0
      real*8 vecPC, vecPM, vecPH
      real*8 trickle, eta
      real*8 beta(nb), ts(nb)
      real*8 beta_cur
      real*8 baseline
      real*8 vecTcalc(nTimes)
      real*8 D_IC, D_IM, D_IH, D_HOS
      real*8 rho                ! reproting rate

      real*8 D_I1M, D_I1C
      real*8 D_I2M, D_I2C

!     The total number of parameters

      integer nparamtot
      
!     Declare items to be returned
      
      real*8 rtn_inf(nTimes,nReals)
      real*8 rtn_pop(nTimes,nReals)
      integer rtn_daily(ndays,nReals)
      integer rec_daily(ndays,nReals)
      integer hos_daily(ndays,nReals)
      integer cli_daily(ndays,nReals)
      integer dead_daily(ndays,nReals)

      integer totalN
      real*8 vecOne
      real*8 tps(ndays+1)
      real*8 epsilon

      real*8 vecHazExI1_M, vecHazExI2_M
      real*8 vecHazExI1_C, vecHazExI2_C
      real*8 vecHazExHOS

      real*8 vecProbExHOS
      real*8 vecProbExI1_C, vecProbExI1_M
      real*8 vecProbExI2_C, vecProbExI2_M
      
      integer I1_M, I1_C, I2_M, I2_C
      integer S, R, HOS, X
      real*8 t_cur, t_next, dt
      integer i, j, k, ind_t
      integer nTimes1, noTPts
      
      real*8 vecFoi,  pVecFoi
      integer noInf, noEntI1_C, noEntI1_M
      integer noEntI2_C, noEntI2_M, noExHOS, noDead
      integer noEntH
      integer noExI1_M, noExI1_C
      integer noExI2_M, noExI2_C 
      integer ignbin
      external ignbin

      ! total number of parameters

      nparamtot = nparam + nb * 2
      
      !     zero all needed arrays

      rtn_inf = 0
      rtn_pop = 0
      rtn_daily = 0
      rec_daily = 0
      hos_daily = 0
      cli_daily = 0
      dead_daily = 0

!     Need to start from day -1 because cases start adding up already on day zero

      nTimes1 = nTimes + 1
      noTPts = nTimes1

      do i = 1, ndays+1
         tps(i) = (i-1)
      enddo

      epsilon = 1.0e-10

      do i = 1, nReals
         
         par = par_tab(i, 1:nparamtot)
!        if (mod(i, 10) == 0)  print*, 'Calculating Realization ', i

!     retrieve parameters from par array


         vecPC = par(1)
         vecPH = par(2)
         D_IC  = par(3)
         D_IM  = par(4)
         D_IH  = par(5)
         D_HOS = par(6)
         rho   = par(7)
         eta   = par(8)
         baseline = par(9)

         do k = 1, nb
            beta(k) = par((k+nparam))
            ts(k)   = par((k+nparam+nb))
         enddo
         
! convert from absolute days in each R(t) value to relative number of days

         do k = 2, (nb-1)
            ts(k) = ts((k-1)) + ts(k)
         enddo
      
!     split each of the two I states to two compartments
!     D_IM/D_IC - the recoverytime from the Im and Ic infectious states - can be the same
!     D_IS - the time it takes to go from clinical to hospital

         D_I1M = D_IM/2.0
         D_I2M = D_IM/2.0
         D_I1C = D_IC/2.0
         D_I2C = D_IC/2.0
         
!     Define housekeeping variables

         totalN = vecN
         dt     = vecTcalc(2) - vecTcalc(1)
         vecOne = 1.0
         vecPM  = vecOne - vecPC

!     Initiate constant hazards
         
         vecHazExI1_M = 1/D_I1M
         vecHazExI2_M = 1/D_I2M
         vecHazExI1_C = 1/D_I1C
         vecHazExI2_C = 1/D_I2C
                  
         vecHazExHOS  = 1/D_HOS !# same rate for exiting for now with just a probability to either go to R/X
!##   Initiate constant probs
         
         vecProbExI1_M = 1 - exp(-dt * vecHazExI1_M)
         vecProbExI2_M = 1 - exp(-dt * vecHazExI2_M)
         vecProbExI1_C = 1 - exp(-dt * vecHazExI1_C)
         vecProbExI2_C = 1 - exp(-dt * vecHazExI2_C)

         vecProbExHOS  = 1 - exp(-dt * vecHazExHOS)

         I1_M = int(vecI0 * vecPM)
         I1_C = int(vecI0 * vecPC)
         
         I2_M = 0
         I2_C = 0
                 
         S = int(vecN * eta - I1_M - I1_C)        
         R = 0
         HOS = 0
         X = 0

         rtn_pop(1,i) = S + I1_M + I2_M + I1_C + I2_C + HOS + R
         
!     initiate the output arrays 
         
!     Initiate the time loop
         j = 2
!     initiate for daily output
         ind_t = 2
         t_next = tps(ind_t)
         
         do while (j <= nTimes)

            t_cur = vecTcalc(j)

            beta_cur = (beta(1) + beta(nb))
            do k = 2, nb
               beta_cur = beta_cur + (beta(k) - beta((k-1))) *
     $              tanh((t_cur-ts((k-1)))/wl)
            enddo
                       
            beta_cur = beta_cur * 0.5
            vecFoi = beta_cur * (I1_M + I2_M + I1_C + I2_C)/vecN +
     $           trickle/(totalN * 7 * dt)

!     Calculate variable probabilites
            pVecFoi = 1 - exp(-dt * vecFoi)
            
            noInf = ignbin(S,  pVecFoi)
            noEntI1_C = ignbin(noInf, vecPC)
            noEntI1_M = ignbin(noInf, vecPM)
            noExI1_M = ignbin(I1_M, vecProbExI1_M)
            noExI2_M = ignbin(I2_M, vecProbExI2_M)
            noExI1_C = ignbin(I1_C, vecProbExI1_C)
            noExI2_C = ignbin(I2_C, vecProbExI2_C)

            noEntH = ignbin(noExI2_C, vecPH)

            noExHOS  = ignbin(HOS, vecProbExHOS)
            

!     # Since we do not keep track of dead/recovered we use the same rate to go
!     # from hospital compartment to recovered or deceased 
!     # we know that for flu the case fatality rate is ~0.1%
            
            noDead = int(noExHOS * 0.03)
!     ## Update the state variables
            
            S = S - noInf
            
            I1_M = I1_M + noEntI1_M - noExI1_M
            I2_M = I2_M + noExI1_M  - noExI2_M
            
            I1_C = I1_C + noEntI1_C - noExI1_C
            I2_C = I2_C + noExI1_C  - noExI2_C
            
            HOS = HOS + noEntH - noExHOS !# noExHos includes recovered and dead 

            R =  R + noExI2_M + noExI2_C + (noExHOS - noDead)
            
!     The dead stage
            X = X + noDead 

!##   Record the other output variables
            rtn_inf(j, i) = noInf
            rtn_pop(j, i) = S + I1_M + I2_M + I1_C + I2_C + HOS + R + X
            
            rtn_daily(ind_t - 1, i) = rtn_daily(ind_t - 1, i) + noInf
            rec_daily(ind_t - 1, i) = rec_daily(ind_t - 1, i) + noExI2_M
     $           + noExI2_C + noExHOS - noDead
            hos_daily(ind_t - 1, i) = hos_daily(ind_t - 1, i) + noEntH
            cli_daily(ind_t - 1, i) = cli_daily(ind_t - 1, i)
     $           + noEntI1_C    !# daily new clinical cases
            dead_daily(ind_t - 1, i) = dead_daily(ind_t - 1, i) + noDead                     

            j = j + 1

            if (t_cur > (t_next + epsilon)) then
               ind_t = ind_t + 1
               t_next = tps(ind_t)
            endif
            
         enddo

         
      enddo

!     multiply hos_daily by the reporting rate

      hos_daily = (hos_daily + baseline) * rho
      
      return

      end subroutine stosirh
      
      
