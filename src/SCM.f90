!
!     SCM, stochastic collocation method
!
!     author  - jure.ravnik@um.si
!
!
!     Include modules
!
PROGRAM StoCMet

      USE GaussHermiteIntegration
      USE GaussLegendreIntegration
      USE SmolyakIntegration

      IMPLICIT NONE

      INTEGER nCP,nRV,nSR,PDFtype,iSR
      REAL, ALLOCATABLE :: RVmu(:),RVsig(:),RVmin(:),RVmax(:)
      CHARACTER(50), ALLOCATABLE :: RVname(:),fmt
      REAL, ALLOCATABLE :: simRuns(:,:),nodes(:,:),simRes(:),simData(:,:)

      INTEGER i,j,ierr
      REAL Ev,Var,skew,kurt



      INTEGER nCol
      CHARACTER(100), ALLOCATABLE :: colNames(:),heName(:)
      REAL, ALLOCATABLE :: EVstatRes(:,:),VARstatRes(:,:),SKEWstatRes(:,:),KURTstatRes(:,:)
      REAL, ALLOCATABLE :: heEV(:),heVAR(:),heSKEW(:),heKURT(:)

      REAL time

      INTEGER nFiles,nSections,iFil,iSec,dataStart,dataLength,fPos,iCol
      CHARACTER(255) fName,path,ime

      CHARACTER(10000) vrstica
      CHARACTER(255),ALLOCATABLE :: bef(:),aft(:)
      INTEGER, PARAMETER :: normal  = 1
      INTEGER, PARAMETER :: uniform = 2
      INTEGER, PARAMETER :: sparseuniform = 3
!
!     Read input file
!
      OPEN (UNIT=12,FILE="scm.inp",STATUS="OLD")

      CALL rOneTL(12,vrstica)
      READ(vrstica,*) nRV
      CALL rOneTL(12,vrstica)
      READ(vrstica,*) nCP
      CALL rOneTL(12,vrstica)
      IF (TRIM(vrstica).EQ."sparseuniform") PDFtype = sparseuniform
      IF (TRIM(vrstica).EQ."uniform") PDFtype = uniform
      IF (TRIM(vrstica).EQ."normal") PDFtype = normal
      IF (TRIM(vrstica).NE."normal".AND.TRIM(vrstica).NE."uniform".AND.TRIM(vrstica).NE."sparseuniform") &
         WRITE(*,*) "ERROR : unknown PDF type!!"
!
!     Names of random variables
!
      ALLOCATE ( RVname(nRV) )

      IF (PDFtype.EQ.uniform.OR.PDFtype.EQ.sparseuniform) THEN

!
!       RV min and max needed
!
        ALLOCATE ( RVmin(nRV) )
        ALLOCATE ( RVmax(nRV) )
!
!       Read
!
        DO i=1,nRV
          CALL rOneTL(12,vrstica)
          READ(vrstica,*) RVname(i),RVmin(i),RVmax(i)
        END DO

      END IF


      IF (PDFtype.EQ.normal) THEN
!
!       RV average value and sigma
!
        ALLOCATE ( RVmu(nRV) )
        ALLOCATE ( RVsig(nRV) )
!
!       Read
!
        DO i=1,nRV
          CALL rOneTL(12,vrstica)
          READ(vrstica,*) RVname(i),RVmu(i),RVsig(i)
        END DO

      END IF
!
!     Init integration routine
!
      IF (PDFtype.EQ.uniform) CALL GLinit(nCP)
      IF (PDFtype.EQ.normal) CALL GHinit(nCP)
      IF (PDFtype.EQ.sparseuniform) CALL Smolyak_init(1,nRV,nCP) ! 5 =Gauss Legendre

!
!     Set up list of parameter values for simulation runs
!
      IF (PDFtype.EQ.sparseuniform) THEN

        nSR = SSQuad%n
        ALLOCATE ( simRuns(nSR,nRV) )
        DO i=1,nSR
          DO j=1,nRV
            simRuns(i,j) =  RVmin(j)+( RVmax(j) - RVmin(j) )*0.5D0*(SSquad%xi(j,i)+1.0D0)  ! uniform distribution of error assumed
          END DO
        END DO

      ELSE

        ALLOCATE (Nodes(nCP,nRV))
        DO i = 1, nCP
          DO j = 1, nRV
            IF (PDFtype.EQ.normal)  Nodes(i,j)=RVmu(j)+GHquad%xi(i)*RVsig(j)*sqrt(2.0D0)  ! natural distribution of error assumed
            IF (PDFtype.EQ.uniform) Nodes(i,j)=RVmin(j)+( RVmax(j) - RVmin(j) )*0.5D0*(GLquad%xi(i)+1.0D0)  ! uniform distribution of error assumed
          END DO
        END DO
!
!       Set up list of parameters for simulation runs
!
        nSR = nCP**nRV
        ALLOCATE ( simRuns(nSR,nRV) )   
        CALL MakeSimRunsList(simRuns,Nodes,nSR,nCP,nRV)

      END IF

!
!     Write simulation runs to file
!
      CALL WriteSimRuns(simRuns,nSR,nRV,RVname)
!
!     Consider analysing the results
!    
      CALL rOneTL(12,vrstica)
      IF (TRIM(vrstica).EQ."brain") THEN
!
!       Analyse results of simulations
!
!
!       Set up array for simulation results
!
        ALLOCATE ( simRes(nSR) )     
!
!       Path to simulation results file folders (1,2,3,... folder names assumed)
!
        CALL rOneTL(12,path)
!
!       Number of files to analyse
!
        CALL rOneTL(12,vrstica)
        READ(vrstica,*) nFiles
!
!       Loop over files to analyse
!
        DO iFil = 1,nFiles
!
!         Name of file to analyse
!
          CALL rOneTL(12,fName)
          WRITE (*,*) TRIM(fName)
!
!         Generate names for result files
!
          WRITE (vrstica,'(A,A)') "ExpV-",TRIM(fname)
          OPEN (UNIT=13,FILE=TRIM(vrstica),STATUS="UNKNOWN")
          WRITE (vrstica,'(A,A)') "Vari-",TRIM(fname)
          OPEN (UNIT=14,FILE=TRIM(vrstica),STATUS="UNKNOWN")
          WRITE (vrstica,'(A,A)') "SDev-",TRIM(fname)
          OPEN (UNIT=15,FILE=TRIM(vrstica),STATUS="UNKNOWN")
          WRITE (vrstica,'(A,A)') "Kurt-",TRIM(fname)
          OPEN (UNIT=16,FILE=TRIM(vrstica),STATUS="UNKNOWN")
          WRITE (vrstica,'(A,A)') "Skew-",TRIM(fname)
          OPEN (UNIT=17,FILE=TRIM(vrstica),STATUS="UNKNOWN")                              
!
!         Number of sections of a file to analyse
!
          CALL rOneTL(12,vrstica)
          READ(vrstica,*) nSections
!
!         Position in file 
!
          fPos = 0
          DO iSec = 1,nSections
!
!           Analyse file
!
            CALL rOneTL(12,vrstica)
            READ(vrstica,*) iCol,dataStart,dataLength
!
!           Copy file parts, which are not analysed to output files
!            
            WRITE (vrstica,'(A,A,A)')  TRIM(path),"1/",TRIM(fname)
            OPEN (UNIT=20,FILE=TRIM(vrstica),STATUS="OLD")
            DO i=1,fPos ! skip to current position
              READ(20,'(A)') vrstica  
            END DO
            DO i=fPos+1,dataStart-1
              READ(20,'(A)') vrstica
              DO j=13,17
                WRITE (j,'(A)') TRIM(vrstica)
              END DO
            END DO
            CLOSE(20)
            fPos = dataStart + dataLength - 1
!  
!           Read from simulation results file
!
            ALLOCATE (simData(nSR,dataLength))
            ALLOCATE (bef(dataLength))
            ALLOCATE (aft(dataLength))

            DO iSR=1,nSR
!  
!             OPEN simulation results file
!
              WRITE (vrstica,'(A,I0,A,A)')  TRIM(path),iSR,"/",TRIM(fname)
              OPEN (UNIT=20,FILE=TRIM(vrstica),STATUS="OLD")
              DO i=1,dataStart-1
                READ(20,'(A)') vrstica  ! dont need this part
              END DO
              DO i=1,dataLength
                READ(20,'(A)') vrstica
                CALL ReadNthCol(vrstica,iCol,simData(iSR,i),bef(i),aft(i))
              END DO
            END DO
!
!           Do statistics
!            
            DO i=1,dataLength
                DO iSR=1,nSR
                  simRes(iSR)=simData(iSR,i)
                END DO
                IF (PDFtype.EQ.normal)  CALL nSCMestStatPar(simRes,nSR,nCP,nRV,Ev,Var,skew,kurt)
                IF (PDFtype.EQ.uniform) CALL uSCMestStatPar(simRes,nSR,nCP,nRV,Ev,Var,skew,kurt)
                IF (PDFtype.EQ.sparseuniform) CALL sparseSCMestStatPar(simRes,nSR,nRV,Ev,Var,skew,kurt)
                WRITE (13,*) TRIM(bef(i)),Ev,TRIM(aft(i))
                WRITE (14,*) TRIM(bef(i)),Var,TRIM(aft(i))
                WRITE (15,*) TRIM(bef(i)),SQRT(Var),TRIM(aft(i))
                WRITE (16,*) TRIM(bef(i)),Kurt,TRIM(aft(i))
                WRITE (17,*) TRIM(bef(i)),Skew,TRIM(aft(i))
            END DO

            DEALLOCATE (simData,bef,aft)

          END DO 
!
!         Copy file footer
!          
          WRITE (vrstica,'(A,A,A)')  TRIM(path),"1/",TRIM(fname)
          OPEN (UNIT=20,FILE=TRIM(vrstica),STATUS="OLD")
          DO i=1,fPos ! skip to current position
            READ(20,'(A)') vrstica  
          END DO
          DO WHILE (.true.)          
            READ(20,'(A)',END=10) vrstica
            DO j=13,17
              WRITE (j,'(A)') TRIM(vrstica)
            END DO
          END DO
10        CONTINUE
          CLOSE(20)
!
!         Close statstics result files
!
          CLOSE(13)
          CLOSE(14)
          CLOSE(15)
          CLOSE(16)                              
          CLOSE(17)          
        END DO ! nFiles

      ELSE IF (TRIM(vrstica).EQ."lio") THEN
!
!       Analyse results of simulations
!
        WRITE (*,*) "Analysing lio results"

!
!       Path to simulation results files
!
        CALL rOneTL(12,path)
        print *,"folder = ",trim(path)
!
!       Get data length        
!
        i=1
        dataLength = 0 
        WRITE (vrstica,'(A,A,I0.4,A)')  TRIM(path),"/rptid",i,".out"
        OPEN (UNIT=20,FILE=TRIM(vrstica),STATUS="OLD")
        CALL rOneTL(20,vrstica)
        DO WHILE (vrstica(1:3).NE."EOF")
          CALL rOneTL(20,vrstica)
          dataLength = dataLength + 1
        END DO
        CLOSE (20)

        print *,"dataLength=",dataLength,nSR
!
!       Allocate space for simulation results   
!                
        ALLOCATE ( simData(nSR,dataLength))  
!
!       Allocate space for statistics results   
!                
        nCol = 15
        ALLOCATE ( colNames(nCol))       
        ALLOCATE ( EVstatRes(dataLength,nCol)) ! Ev,Var,skew,kurt
        ALLOCATE ( VARstatRes(dataLength,nCol)) ! Ev,Var,skew,kurt
        ALLOCATE ( SKEWstatRes(dataLength,nCol)) ! Ev,Var,skew,kurt
        ALLOCATE ( KURTstatRes(dataLength,nCol)) ! Ev,Var,skew,kurt
!
!       Time 
!
        iCol = 1
        colNames(icol)="Time_[min]" ! time is not analysed
        write (*,*) "Working on: ",trim(colNames(icol))
        WRITE (vrstica,'(A,A,I0.4,A)')  TRIM(path),"/rptid",i,".out"
        OPEN (UNIT=20,FILE=TRIM(vrstica),STATUS="OLD")
        DO i=1,dataLength
          CALL rOneTL(20,vrstica)
          READ(vrstica,*) time
          EVstatRes(i,iCol) = time
          VARstatRes(i,iCol) = time
          SKEWstatRes(i,iCol) = time
          KURTstatRes(i,iCol) = time
        END DO
        CLOSE (20)
!
!       file = rptid, column = 2
!
        iCol = iCol + 1
        colNames(iCol)="Rp_[cm2_torr_h/g]"
        write (*,*) "Working on: ",trim(colNames(icol))
        call readSimData(nSR,path,"/rptid",2,simData,dataLength)
        call doStatistics(simData,dataLength,nSR,nCP,nRV,EVstatRes,VARstatRes,SKEWstatRes,KURTstatRes,nCol,iCol,PDFtype)
!
!       file = jurid, column = 2
!
        iCol = iCol + 1
        colNames(iCol)="T[C]_at_y=0mm"
        write (*,*) "Working on: ",trim(colNames(icol))
        call readSimData(nSR,path,"/jurid",2,simData,dataLength)
        call doStatistics(simData,dataLength,nSR,nCP,nRV,EVstatRes,VARstatRes,SKEWstatRes,KURTstatRes,nCol,iCol,PDFtype)
!
!       file = jurid, column = 3
!
        iCol = iCol + 1
        colNames(iCol)="T[C]_at_y=0.5mm"
        write (*,*) "Working on: ",trim(colNames(icol))
        call readSimData(nSR,path,"/jurid",3,simData,dataLength)
        call doStatistics(simData,dataLength,nSR,nCP,nRV,EVstatRes,VARstatRes,SKEWstatRes,KURTstatRes,nCol,iCol,PDFtype)
!
!       file = jurid, column = 4
!
        iCol = iCol + 1
        colNames(iCol)="T[C]_at_y=1.0mm"
        write (*,*) "Working on: ",trim(colNames(icol))
        call readSimData(nSR,path,"/jurid",4,simData,dataLength)
        call doStatistics(simData,dataLength,nSR,nCP,nRV,EVstatRes,VARstatRes,SKEWstatRes,KURTstatRes,nCol,iCol,PDFtype)
!
!       file = jurid, column = 5
!
        iCol = iCol + 1
        colNames(iCol)="T[C]_at_y=1.5mm"
        write (*,*) "Working on: ",trim(colNames(icol))
        call readSimData(nSR,path,"/jurid",5,simData,dataLength)
        call doStatistics(simData,dataLength,nSR,nCP,nRV,EVstatRes,VARstatRes,SKEWstatRes,KURTstatRes,nCol,iCol,PDFtype)
!
!       file = jurid, column = 6
!
        iCol = iCol + 1
        colNames(iCol)="T[C]_at_y=2.0mm"
        write (*,*) "Working on: ",trim(colNames(icol))
        call readSimData(nSR,path,"/jurid",6,simData,dataLength)
        call doStatistics(simData,dataLength,nSR,nCP,nRV,EVstatRes,VARstatRes,SKEWstatRes,KURTstatRes,nCol,iCol,PDFtype)
!
!       file = jurid, column = 7
!
        iCol = iCol + 1
        colNames(iCol)="T[C]_at_y=2.5mm"
        write (*,*) "Working on: ",trim(colNames(icol))
        call readSimData(nSR,path,"/jurid",7,simData,dataLength)
        call doStatistics(simData,dataLength,nSR,nCP,nRV,EVstatRes,VARstatRes,SKEWstatRes,KURTstatRes,nCol,iCol,PDFtype)
!
!       file = jurid, column = 8
!
        iCol = iCol + 1
        colNames(iCol)="Tpar"
        write (*,*) "Working on: ",trim(colNames(icol))
        call readSimData(nSR,path,"/jurid",8,simData,dataLength)
        call doStatistics(simData,dataLength,nSR,nCP,nRV,EVstatRes,VARstatRes,SKEWstatRes,KURTstatRes,nCol,iCol,PDFtype)
!
!       file = jurid, column = 9
!
        iCol = iCol + 1
        colNames(iCol)="dmdt[g/h]"
        write (*,*) "Working on: ",trim(colNames(icol))
        call readSimData(nSR,path,"/jurid",9,simData,dataLength)
        call doStatistics(simData,dataLength,nSR,nCP,nRV,EVstatRes,VARstatRes,SKEWstatRes,KURTstatRes,nCol,iCol,PDFtype)
!
!       file = jurid, column = 10
!
        iCol = iCol + 1
        colNames(iCol)="%dried"
        write (*,*) "Working on: ",trim(colNames(icol))
        call readSimData(nSR,path,"/jurid",10,simData,dataLength)
        call doStatistics(simData,dataLength,nSR,nCP,nRV,EVstatRes,VARstatRes,SKEWstatRes,KURTstatRes,nCol,iCol,PDFtype)
!
!       file = jurid, column = 11
!
        iCol = iCol + 1
        colNames(iCol)="RpBEM"
        write (*,*) "Working on: ",trim(colNames(icol))
        call readSimData(nSR,path,"/jurid",11,simData,dataLength)
        call doStatistics(simData,dataLength,nSR,nCP,nRV,EVstatRes,VARstatRes,SKEWstatRes,KURTstatRes,nCol,iCol,PDFtype)
!
!       file = jurid, column = 12
!
        iCol = iCol + 1
        colNames(iCol)="Tmax"
        write (*,*) "Working on: ",trim(colNames(icol))
        call readSimData(nSR,path,"/jurid",12,simData,dataLength)
        call doStatistics(simData,dataLength,nSR,nCP,nRV,EVstatRes,VARstatRes,SKEWstatRes,KURTstatRes,nCol,iCol,PDFtype)
!
!       file = jurid, column = 13
!
        iCol = iCol + 1
        colNames(iCol)="T[C]_at_y=5.5mm"
        write (*,*) "Working on: ",trim(colNames(icol))
        call readSimData(nSR,path,"/jurid",13,simData,dataLength)
        call doStatistics(simData,dataLength,nSR,nCP,nRV,EVstatRes,VARstatRes,SKEWstatRes,KURTstatRes,nCol,iCol,PDFtype)
!
!       file = jurid, column = 14
!
        iCol = iCol + 1
        colNames(iCol)="T[C]_at_y=10.5mm"
        write (*,*) "Working on: ",trim(colNames(icol))
        call readSimData(nSR,path,"/jurid",14,simData,dataLength)
        call doStatistics(simData,dataLength,nSR,nCP,nRV,EVstatRes,VARstatRes,SKEWstatRes,KURTstatRes,nCol,iCol,PDFtype)                
!
!       PD time
!
        ALLOCATE ( simRes(nSR) )
        simRes = 0.0D0
        DO iSR=1,nSR ! ,5 BRISI
          WRITE (vrstica,'(A,A,I0.4,A)')  TRIM(path),"/jurid",iSR,".out"          
          OPEN (UNIT=20,FILE=TRIM(vrstica),STATUS="OLD")
          DO i=1,dataLength
            CALL rOneTL(20,vrstica)
          END DO
          CALL rOneTL(20,vrstica)
          read(vrstica,*) simRes(iSR),simRes(iSR)
          CLOSE (20)         
        END DO
        IF (PDFtype.EQ.normal) CALL nSCMestStatPar(simRes,nSR,nCP,nRV,Ev,Var,skew,kurt)
        IF (PDFtype.EQ.uniform) CALL uSCMestStatPar(simRes,nSR,nCP,nRV,Ev,Var,skew,kurt)
        IF (PDFtype.EQ.sparseuniform) CALL sparseSCMestStatPar(simRes,nSR,nRV,Ev,Var,skew,kurt)
     
        deallocate(simRes)        
!
!       Generate names for result files
!

        OPEN (UNIT=13,FILE=TRIM("scm-results.dat"),STATUS="UNKNOWN")
!       title line
        WRITE(13,"(A)") "# SCM result file"
        fmt=""
        write(fmt,*) Ev
        WRITE(13,"(A)") "# PDTIME ExpV = "//fmt
        write(fmt,*) Var
        WRITE(13,"(A)") "# PDTIME Vari = "//fmt
        write(fmt,*) sqrt(Var)
        WRITE(13,"(A)") "# PDTIME SDev = "//fmt
        write(fmt,*) Skew
        WRITE(13,"(A)") "# PDTIME Skew = "//fmt
        write(fmt,*) Kurt
        WRITE(13,"(A)") "# PDTIME Kurt = "//fmt


        i=1
        write(13,"(A,I0,1X,A)") "# ",i,TRIM(colNames(1))
        do j=2,nCol
          i=i+1
          write(13,"(A,I0,1X,A)") "# col ",i," = ExpV_"//TRIM(colNames(j))
          i=i+1
          write(13,"(A,I0,1X,A)") "# col ",i," = Vari_"//TRIM(colNames(j))
          i=i+1
          write(13,"(A,I0,1X,A)") "# col ",i," = SDev_"//TRIM(colNames(j))
          i=i+1
          write(13,"(A,I0,1X,A)") "# col ",i," = Skew_"//TRIM(colNames(j))                
          i=i+1
          write(13,"(A,I0,1X,A)") "# col ",i," = Kurt_"//TRIM(colNames(j))                        
        end do
        
        write(fmt,'(A1,I0,A2)') "(",5*(nCol-1)+2,"A)"
        write(13,fmt) "# ",TRIM(colNames(1)),( &
          " ExpV_"//TRIM(colNames(i))// &
          " Vari_"//TRIM(colNames(i))// & 
          " SDev_"//TRIM(colNames(i))// &
          " Skew_"//TRIM(colNames(i))// &
          " Kurt_"//TRIM(colNames(i)) &
          ,i=2,nCol)

        do i=1,dataLength
          write(13,*)  EVstatRes(i,1), &  ! time
                       (EVstatRes(i,j),VARstatRes(i,j),sqrt(VARstatRes(i,j)),SKEWstatRes(i,j),KURTstatRes(i,j),j=2,nCol)
        end do
!
!       Close statstics result files
!
        CLOSE(13)

      ELSE IF (TRIM(vrstica).EQ."heatEx") THEN
        WRITE (*,*) "Heat exchanger postprocessing, ...."

!
!       Path to simulation results files
!
        CALL rOneTL(12,path)
        print *,"folder = ",trim(path)
!
!       Get data length        
!

!       data lenght = koliko razlicnih "rezultatov" analiziramo, nSR = število random spremenljivk        
        WRITE (fName,'(A,A)')  TRIM(path),"/Results.dat"
        OPEN (UNIT=20,FILE=TRIM(fName),STATUS="OLD")
        CALL rOneTL(20,vrstica)
        CALL countCols( vrstica, dataLength )
        CLOSE (20)

        print *,"dataLength=",dataLength,nSR
!
!       Allocate space for simulation results   
!                
        ALLOCATE ( simData(nSR,dataLength))  

        OPEN (UNIT=20,FILE=TRIM(fName),STATUS="OLD")
        DO i = 1,nSR
                CALL rOneTL(20,vrstica)
                DO j=1,dataLength
                        CALL ReadNthColSimple(vrstica,j,simData(i,j))
                END DO
        END DO
        CLOSE(20)

        ALLOCATE ( simRes(nSR) )
        ALLOCATE (heEV(dataLength))
        ALLOCATE (heVAR(dataLength))
        ALLOCATE (heSKEW(dataLength))
        ALLOCATE (heKURT(dataLength))
        ALLOCATE (heName(dataLength))

        OPEN (UNIT=20,FILE=TRIM(fName),STATUS="OLD")
        READ(20,*) heName
        CLOSE(20)

        OPEN (UNIT=13,FILE=TRIM("scm-results.dat"),STATUS="UNKNOWN")
!       title line
        WRITE(13,"(A)") "# SCM result file"

        WRITE(13,"(A)") "# Name expected_value st_dev variance skewness kurtosis"

        DO i=1,dataLength
          DO iSR=1,nSR
            simRes(iSR)=simData(iSR,i)
          END DO
          IF (PDFtype.EQ.normal) CALL             nSCMestStatPar(simRes,nSR,nCP,nRV,heEV(i),heVAR(i),heSKEW(i),heKURT(i))
          IF (PDFtype.EQ.uniform) CALL            uSCMestStatPar(simRes,nSR,nCP,nRV,heEV(i),heVAR(i),heSKEW(i),heKURT(i))
          IF (PDFtype.EQ.sparseuniform) CALL sparseSCMestStatPar(simRes,nSR,nRV,    heEV(i),heVAR(i),heSKEW(i),heKURT(i))
        
          print *,i-1,trim(heName(i))
          WRITE (ime,'(A,I0,A,A)') "_",i-1,"_",trim(heName(i))
          WRITE (13,*) trim(ime),heEV(i),sqrt(heVAR(i)),heVAR(i),heSKEW(i),heKURT(i)
        
        END DO
      
        deallocate(simRes,heEV,heVAR,heSKEW,heKURT,heName)
        CLOSE(13)

      ELSE IF (TRIM(vrstica).EQ."bioheat") THEN
!
!       Analyse results of simulations
!
        WRITE (*,*) "Analysing bioheat OpenFOAM results"

!
!       Path to simulation results files
!
        CALL rOneTL(12,path)
        print *,"folder = ",trim(path)
        print *,"nSR = ",nSR
!
!       Set up array for simulation results
!
        dataLength = 10  ! 2141995 ! = number of lines in result files all_T_values.txt  
        ALLOCATE ( simRes(nSR) )      
        ALLOCATE (heEV(dataLength))  ! Expected value
        ALLOCATE (heVAR(dataLength)) ! Variance
        ALLOCATE (heSKEW(dataLength)) ! Skewness
        ALLOCATE (heKURT(dataLength)) ! Kurtosis
!
!       Set up results file
!        
        OPEN (UNIT=13,FILE=TRIM("scm-results-EV.dat"),STATUS="UNKNOWN")
        OPEN (UNIT=14,FILE=TRIM("scm-results-STDEV.dat"),STATUS="UNKNOWN")
        OPEN (UNIT=15,FILE=TRIM("scm-results-VAR.dat"),STATUS="UNKNOWN")
        OPEN (UNIT=16,FILE=TRIM("scm-results-SKEW.dat"),STATUS="UNKNOWN")
        OPEN (UNIT=17,FILE=TRIM("scm-results-KURT.dat"),STATUS="UNKNOWN")
!
!       Loop over lines in files
!                    
        DO i = 1, dataLength ! number of lines in result files
!           Loop over files to analyse                
            DO iFil = 1,nSR
                WRITE(fname,'(A,A,I0,A)') TRIM(path),"/",iFil,"/all_T_values.txt"
                CALL ReadNumberFromLine(fname,i,simRes(iFil),ierr)
                if (ierr.NE.0) THEN
                        WRITE(*,*) "Error reading file: ",TRIM(fname)," line: ",i
                        STOP
                END IF
            END DO
!           Do statistics              
            IF (PDFtype.EQ.normal) CALL             nSCMestStatPar(simRes,nSR,nCP,nRV,heEV(i),heVAR(i),heSKEW(i),heKURT(i))
            IF (PDFtype.EQ.uniform) CALL            uSCMestStatPar(simRes,nSR,nCP,nRV,heEV(i),heVAR(i),heSKEW(i),heKURT(i))
            IF (PDFtype.EQ.sparseuniform) CALL sparseSCMestStatPar(simRes,nSR,nRV,    heEV(i),heVAR(i),heSKEW(i),heKURT(i))  
!           Save results to file
            WRITE (13,*) heEV(i)
            WRITE (14,*) sqrt(heVAR(i))
            WRITE (15,*) heVAR(i)
            WRITE (16,*) heSKEW(i)
            WRITE (17,*) heKURT(i)
        END DO

        deallocate(simRes,heEV,heVAR,heSKEW,heKURT)
        CLOSE(13)
        CLOSE(14)
        CLOSE(15)
        CLOSE(16)
        CLOSE(17)
        
      ELSE 
        WRITE (*,*) "No postprocessing, quitting!"
      END IF

!
!     Close input file
!
      CLOSE (12)


      END

