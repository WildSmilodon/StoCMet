
      SUBROUTINE ReadNumberFromLine(fileName,n,value,ierr)

      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: fileName

      INTEGER, INTENT(IN) :: n
      REAL, INTENT(OUT) :: value
      INTEGER, INTENT(OUT) :: ierr

      INTEGER :: i,ios
      INTEGER, PARAMETER :: unitNo = 21

      value = 0.0D0
      ierr = 0

      IF (n.LT.1) THEN
        ierr = 1
        RETURN
      END IF

      OPEN(UNIT=unitNo,FILE=TRIM(fileName),STATUS="OLD",ACTION="READ",IOSTAT=ios)
      IF (ios.NE.0) THEN
        ierr = 2
        RETURN
      END IF

      DO i=1,n
        READ(unitNo,*,IOSTAT=ios) value
        IF (ios.NE.0) THEN
          ierr = 3
          CLOSE(unitNo)
          RETURN
        END IF
      END DO

      CLOSE(unitNo,IOSTAT=ios)
      IF (ios.NE.0) ierr = 4

      END SUBROUTINE ReadNumberFromLine




!
!
!
!
!  ------------------------------------------------------------------
!
      SUBROUTINE MakeSimRunsList(simRuns,Nodes,nSR,nCP,nRV)
   
      IMPLICIT NONE

      INTEGER, ALLOCATABLE :: ind(:)
      REAL simRuns(nSR,nRV)
      REAL nodes(nCP,nRV)
      INTEGER i,j,ii,nSR,nCP,nRV
      LOGICAL done      

      ALLOCATE ( ind(nRV) )
      ind = 1 
      DO i = 1,nSR
         DO j=1,nRV
           simRuns(i,j) = Nodes(ind(j),j)
         END DO

         done = .FALSE.
         ii = nRV
         DO WHILE (.not.done.AND.i.NE.nSR)
  
           IF (ind(ii).LT.nCP) THEN
             ind(ii) =  ind(ii) + 1
             done = .TRUE.
           ELSE
             ind(ii) = 1
             ii = ii - 1
           END IF
         END DO
      END DO

      END SUBROUTINE


!
!  ------------------------------------------------------------------
!
      SUBROUTINE WriteSimRuns(simRuns,nSR,nRV,RVname)
   
      IMPLICIT NONE

      INTEGER i,j,nSR,nRV
      REAL simRuns(nSR,nRV)
      CHARACTER*(*) RVname(nRV)

!      WRITE (*,*) (TRIM(RVname(j))//" ",j=1,nRV)
!      DO i = 1,nSR
!         WRITE (*,*) (simRuns(i,j),j=1,nRV)
!      END DO

      WRITE (*,*) "Number of random variables :",nRV
      WRITE (*,*) "Number of runs of deterministic code :",nSR

      OPEN (UNIT=53,FILE=TRIM("scm-runList.txt"),STATUS="UNKNOWN")
      WRITE (53,'(A,I0,A,I0)') " # number of random variables = ",nRV,", number of runs = ",nSR
      WRITE (53,*) "# ID ",(TRIM(RVname(j))//" ",j=1,nRV)
      DO i = 1,nSR
         WRITE (53,*) i,(simRuns(i,j),j=1,nRV)
      END DO

      CLOSE (53)

      END SUBROUTINE


!
!  ------------------------------------------------------------------
!
      SUBROUTINE nSCMestStatPar(simRes,nSR,nCP,nRV,Ev,Var,skew,kurt)
!
!     Estimates Ev, variance, skewness and kurtisis 
!
      USE GaussHermiteIntegration      
      IMPLICIT NONE
      INTEGER nSR,nRV,i,j,nCP,ii
      INTEGER, ALLOCATABLE :: ind(:)
      REAL pi,Ev,Var,simRes(nSR),w,s3,s4,s2,s1,skew,kurt
      LOGICAL done

      pi=ATAN(1.0)*4.0D0

      s1 = 0.0D0
      s3 = 0.0D0
      s4 = 0.0D0
      s2 = 0.0D0
      ALLOCATE ( ind(nRV) )  ! tu so indeksi za zanke ( nRV nested loops from 1 to nCP )
      
      ind = 1 
      DO i = 1,nSR
         w = 1.0D0
         DO j=1,nRV
           w = w * GHquad%wi(ind(j))/sqrt(pi)
         END DO

         s1 = s1 + w * simRes(i)
         s2 = s2 + w * simRes(i) * simRes(i)
         s3 = s3 + w * simRes(i) * simRes(i) * simRes(i)
         s4 = s4 + w * simRes(i) * simRes(i) * simRes(i) * simRes(i)

         done = .FALSE.
         ii = nRV

!        ( nRV nested loops from 1 to nCP )         
         DO WHILE (.not.done.AND.i.NE.nSR)  
           IF (ind(ii).LT.nCP) THEN
             ind(ii) =  ind(ii) + 1
             done = .TRUE.
           ELSE
             ind(ii) = 1
             ii = ii - 1
           END IF
         END DO
      END DO

      Ev = s1
      Var = s2 - Ev * Ev
      IF (Var.GT.0.0D0) THEN
        skew = ( s3 - 3.0D0 * Ev * Var - Ev*Ev*Ev)/(Var*Sqrt(Var))
        kurt = ( s4 - 4.0D0 * Ev * s3 + 6.0D0 * Ev * Ev * s2 - 3.0D0 * Ev * Ev * Ev * Ev) / (Var*Var)
      ELSE
        skew = 0.0D0
        kurt = 0.0D0
      END IF


      DEALLOCATE (ind)


      END SUBROUTINE

!
!  ------------------------------------------------------------------
!
      SUBROUTINE uSCMestStatPar(simRes,nSR,nCP,nRV,Ev,Var,skew,kurt)
!
!     Estimates Ev, variance, skewness and kurtisis 
!
      USE GaussLegendreIntegration      
      IMPLICIT NONE
      INTEGER nSR,nRV,i,j,nCP,ii
      INTEGER, ALLOCATABLE :: ind(:)
      REAL pi,Ev,Var,simRes(nSR),w,s3,s4,s2,s1,skew,kurt
      LOGICAL done

      pi=ATAN(1.0)*4.0D0

      s1 = 0.0D0
      s3 = 0.0D0
      s4 = 0.0D0
      s2 = 0.0D0
      ALLOCATE ( ind(nRV) )  ! tu so indeksi za zanke ( nRV nested loops from 1 to nCP )
    
      ind = 1 
      DO i = 1,nSR
         w = 1.0D0
         DO j=1,nRV
           w = w * GLquad%wi(ind(j)) * 0.5D0 
         END DO

         
         s1 = s1 + w * simRes(i)
         s2 = s2 + w * simRes(i) * simRes(i)
         s3 = s3 + w * simRes(i) * simRes(i) * simRes(i)
         s4 = s4 + w * simRes(i) * simRes(i) * simRes(i) * simRes(i)

         done = .FALSE.
         ii = nRV

!        ( nRV nested loops from 1 to nCP )         
         DO WHILE (.not.done.AND.i.NE.nSR)  
           IF (ind(ii).LT.nCP) THEN
             ind(ii) =  ind(ii) + 1
             done = .TRUE.
           ELSE
             ind(ii) = 1
             ii = ii - 1
           END IF
         END DO
      END DO

      Ev = s1
      Var = s2 - Ev * Ev      

      IF (Var.GT.0.0D0) THEN
        skew = ( s3 - 3.0D0 * Ev * Var - Ev*Ev*Ev)/(Var*Sqrt(Var))
        kurt = ( s4 - 4.0D0 * Ev * s3 + 6.0D0 * Ev * Ev * s2 - 3.0D0 * Ev * Ev * Ev * Ev) / (Var*Var)
      ELSE
        Var = 0.0D0
        skew = 0.0D0
        kurt = 0.0D0
      END IF


      DEALLOCATE (ind)


      END SUBROUTINE


!
!  ------------------------------------------------------------------
!
      SUBROUTINE sparseSCMestStatPar(simRes,nSR,nRV,Ev,Var,skew,kurt)
!
!     Estimates Ev, variance, skewness and kurtosis 
!
      USE SmolyakIntegration      
      IMPLICIT NONE
      INTEGER nSR,i,j,nRV
      REAL Ev,Var,simRes(nSR),w,s3,s4,s2,s1,skew,kurt,ep

      s1 = 0.0D0
      s2 = 0.0D0
      s3 = 0.0D0
      s4 = 0.0D0

      ep = 1.0D0 ! (1/2)^nRV
      DO j=1,nRV
        ep = ep * 0.5D0 
      END DO
   
      DO i = 1,nSR
         w = ep * SSquad%wi(i)        
         s1 = s1 + w * simRes(i)
         s2 = s2 + w * simRes(i) * simRes(i)
         s3 = s3 + w * simRes(i) * simRes(i) * simRes(i)
         s4 = s4 + w * simRes(i) * simRes(i) * simRes(i) * simRes(i)

      END DO

      Ev = s1
      Var = s2 - Ev * Ev      

      IF (Var.LT.0.0D0) THEN  ! na novo izracunaj na drugi nacin (morda ne bo negativno, ceprav ponavadi je)
        Var = 0.0D0
        DO i = 1,nSR
           w = ep * SSquad%wi(i)        
           Var = Var + w * ( simRes(i) - Ev ) * ( simRes(i) - Ev )
        END DO
      END IF



      IF (Var.GT.0.0D0) THEN
        skew = ( s3 - 3.0D0 * Ev * Var - Ev*Ev*Ev)/(Var*Sqrt(Var))
        kurt = ( s4 - 4.0D0 * Ev * s3 + 6.0D0 * Ev * Ev * s2 - 3.0D0 * Ev * Ev * Ev * Ev) / (Var*Var)
      ELSE
        Var = 0.0D0    ! varianca vcasih -1.0E-10  zaradi odstevanja podobnih stevil, glej : https://en.wikipedia.org/wiki/Variance
                       ! oziroma zaradi slabe integracije.
        skew = 0.0D0
        kurt = 0.0D0
      END IF

      END SUBROUTINE


!
!  ------------------------------------------------------------------
!
      SUBROUTINE ReadNthCol(line,col,res,bef,aft)
      IMPLICIT NONE
      REAL res
      CHARACTER*(*) line,bef,aft
      INTEGER col,ncol,sIDX,eIDX
      REAL,ALLOCATABLE :: tmp(:)
      CHARACTER(255), ALLOCATABLE :: stmp(:)

      CALL GetNCol(line,ncol)

      ALLOCATE (tmp(nCol),stmp(ncol))

      READ(line,*) tmp
      READ(line,*) stmp
      res=tmp(col)

      IF (ncol.EQ.1) THEN
        bef=" "
        aft = " "
      ELSE IF (col.EQ.1) THEN
        eIDX = INDEX(TRIM(line),TRIM(stmp(col+1)))
        bef=" "
        aft = " "//line(eIDX:len_trim(line))
      ELSE IF (col.LT.ncol) THEN
        sIDX = INDEX(TRIM(line),TRIM(stmp(col)))
        eIDX = INDEX(TRIM(line),TRIM(stmp(col+1)))
        bef = line(1:sIDX-1)//" "
        aft = " "//line(eIDX:len_trim(line))
      ELSE 
        sIDX = INDEX(TRIM(line),TRIM(stmp(col)))
        bef = line(1:sIDX-1)//" "
        aft = " "
      END IF

      DEALLOCATE(tmp)

      END

!
!  ------------------------------------------------------------------
!
subroutine readSimData(nSR,path,name,col,simData,dataLength)
  character*(*) path,name
  character(500) vrstica
  integer nSR,col,dataLength
  real simData(nSR,dataLength)
  integer iSr

  simData = 0.0D0

  do iSR=1,nSR ! ,5 BRISI
!  
! OPEN simulation results file
!
    WRITE (vrstica,'(A,A,I0.4,A)')  TRIM(path),trim(name),iSR,".out"          
!    print *,TRIM(vrstica)
    OPEN (UNIT=20,FILE=TRIM(vrstica),STATUS="OLD",ERR = 100)
    DO i=1,dataLength
      CALL rOneTL(20,vrstica)
      CALL ReadNthColSimple(vrstica,col,simData(iSR,i))       
    END DO
    CLOSE (20)
  END DO   

  return

100   CONTINUE
  print *,"error opening file!"   
  stop
end subroutine    
!
!  ------------------------------------------------------------------
!
subroutine doStatistics(simData,dataLength,nSR,nCP,nRV,EV,VAR,SKEW,KURT,nCol,col,PDFtype)
  implicit none
  integer nSR,col,dataLength,nCol,nCP,nRV
  real simData(nSR,dataLength)
  real EV(dataLength,nCol)
  real VAR(dataLength,nCol)
  real SKEW(dataLength,nCol)
  real KURT(dataLength,nCol)
  real, allocatable :: simRes(:)  
  integer i,iSR,PDFtype
  INTEGER, PARAMETER :: normal  = 1
  INTEGER, PARAMETER :: uniform = 2
  INTEGER, PARAMETER :: sparseuniform = 3  
!
! Set up array for simulation results
!
  ALLOCATE ( simRes(nSR) )

  DO i=1,dataLength
    DO iSR=1,nSR
      simRes(iSR)=simData(iSR,i)
    END DO
    IF (PDFtype.EQ.normal) CALL nSCMestStatPar(simRes,nSR,nCP,nRV,EV(i,col),VAR(i,col),SKEW(i,col),KURT(i,col))
    IF (PDFtype.EQ.uniform) CALL uSCMestStatPar(simRes,nSR,nCP,nRV,EV(i,col),VAR(i,col),SKEW(i,col),KURT(i,col))
    IF (PDFtype.EQ.sparseuniform) CALL sparseSCMestStatPar(simRes,nSR,nRV,EV(i,col),VAR(i,col),SKEW(i,col),KURT(i,col))
  END DO

  deallocate(simRes)

end subroutine


!
!  ------------------------------------------------------------------
!
      SUBROUTINE ReadNthColSimple(line,col,res)
        IMPLICIT NONE
        REAL res
        CHARACTER*(*) line
        INTEGER col,ncol
        REAL,ALLOCATABLE :: tmp(:)
  
        CALL GetNCol(line,ncol)
 
        ALLOCATE (tmp(nCol))
  
        READ(line,*) tmp
        res=tmp(col)
  
        DEALLOCATE(tmp)
  
        END
  
  

!
!  ------------------------------------------------------------------
!
      SUBROUTINE GetNCol(line,ncol)
      IMPLICIT NONE
      CHARACTER*(*) line
      INTEGER ncol,i
      LOGICAL lastSpace

      ncol = 0

      lastSpace = .TRUE.
      DO i=1,LEN_TRIM(line)        
        IF ((line(i:i).EQ." ")) THEN !.AND.(lastSpace.EQV..FALSE.)) THEN 
          lastSpace = .TRUE.
        END IF

        IF (line(i:i).NE." ".AND.lastSpace) THEN
           nCol = nCol + 1
           lastSpace = .FALSE.
        END IF

      END DO

      END


!
!  ------------------------------------------------------------------
!
      SUBROUTINE simulation(x,res)
      IMPLICIT NONE
      REAL x,res

      res = 2.0D0 * x * x  + 3.0D0

      END SUBROUTINE

!
!  ------------------------------------------------------------------
!
      SUBROUTINE simulation2(x1,x2,res)
      IMPLICIT NONE
      REAL x1,x2,res

      res = 2.0D0 * x1 * x1  + 3.0D0 + x2*x2

      END SUBROUTINE

!
!  ------------------------------------------------------------------
!
      SUBROUTINE simulationN(n,x,nX,res)
      IMPLICIT NONE
      INTEGER nX,n
      REAL x(nX),res

      IF (n.EQ.1) THEN
        res = 2.0D0 * x(1) * x(1)  + 3.0D0
      END IF
      IF (n.EQ.2) THEN
        res = 2.0D0 * x(1) * x(1)  + 3.0D0 + x(2)*x(2)  
      END IF
      IF (n.EQ.3) THEN
        res = 2.0D0 * x(1) * x(1)  + 3.0D0 + x(2)*x(2)  + x(3)*x(3)
      END IF

      END SUBROUTINE


!
!  ------------------------------------------------------------------
!


      SUBROUTINE test()

      USE GaussHermiteIntegration

      IMPLICIT NONE

      REAL, ALLOCATABLE :: LagNodes(:) 
      REAL mu,sig,vsota,res,pi
      REAL Ev,var
      INTEGER i

      mu=1.0D0
      sig=0.1D0
      pi=ATAN(1.0)*4.0D0
      print *,pi

       ALLOCATE (LagNodes(GHquad%n))
        DO i = 1, GHquad%n
            LagNodes(i)=mu+GHquad%xi(i)*sig*sqrt(2.0D0)
        END DO

        print *,"GHquad%n",GHquad%n

        print *,"LagNodes",LagNodes

        vsota = 0.0D0
        DO i = 1,GHquad%n 
              CALL simulation(LagNodes(i),res)
              vsota = vsota + res*GHquad%wi(i)
        END DO
        Ev=vsota/sqrt(pi)

        vsota = 0.0D0
        DO i = 1,GHquad%n 
              CALL simulation(LagNodes(i),res)
              vsota = vsota + res*res*GHquad%wi(i)
        END DO
        var=vsota/sqrt(pi) - Ev*Ev


        print *,Ev,var

      END SUBROUTINE

!
!  ------------------------------------------------------------------
!
      SUBROUTINE LagrangeIF(x,xi,n)

      IMPLICIT NONE

      REAL x ! tocka v kateri racunam Lagrangeovo funkcijo
      INTEGER n ! stevilo Lagrangevi tock
      REAL xi(n) ! Lagrangeove tocke
      REAL LagF(n) ! REzultat, vrednosti Langrageoveih interpolacijskih funkcij

      INTEGER i,j

      DO i = 1,n 
        LagF(i)=1.0D0
        DO j = 1,n
          IF (i.NE.j) THEN
            LagF(i) =  LagF(i) * ( x - xi(j) ) / ( xi(i) - xi(j) ) 
          END IF
        END DO
        print *, i,LagF(i)
      END DO

      END SUBROUTINE


!______________________________________________________________________C
!______________________________________________________________________C
      SUBROUTINE rOneTL(lun,OneLine)
!     _    ___ _    _
!     Read One Text Line
!
!______________________________________________________________________C
!     Returns the first nonempty text line in file LUN, which does not
!     include the # character. If end of file is encoutered, it returns EOF
      CHARACTER*(*) OneLine
      INTEGER lun,i

10    READ(lun,'(A)',END=20) OneLine

!     Check if line is empty
      IF (len_trim(OneLine).EQ.0) GOTO 10

!     Check if line contains # character
      DO i=1,len_trim(OneLine)
        IF (OneLine(i:i).EQ.'#') GOTO 10
      ENDDO

      RETURN

20    OneLine='EOF'
      END
      

     subroutine countCols( line, n )
        implicit none
        character(*), intent(in) :: line
        real*8  :: buf( 10000 )
        integer :: n
    
        n = 1
        do
            read( line, *, end=100, err=100 ) buf( 1 : n )   !! (See Appendix for why buf is used here)
            n = n + 1
        enddo
    100 continue
        n = n - 1
    end

