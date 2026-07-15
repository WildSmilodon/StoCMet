C
C     Gauss-Hermite integration library
C     (c) Jure Ravnik, 2019
C
C
C     used for integration of f(x)*exp(-x^2) from -inf to inf
C
C     usage : integral = sum_i ( f(xi)*wi )
C
      MODULE GaussHermiteIntegration
      IMPLICIT NONE
C
C -----------------------------------------------------------------------------------------
C
C
C
C -----------------------------------------------------------------------------------------
C
      TYPE GaussHermiteType
         INTEGER :: n    
         REAL, POINTER :: xi(:),wi(:)
      END TYPE
C
C     Variables
C
      TYPE(GaussHermiteType) :: GHquad
C
C     Subroutines
C
      CONTAINS

C
C -----------------------------------------------------------------------------------------
C
      SUBROUTINE GHtest

      IMPLICIT NONE

      REAL anal
      REAL vsota
      INTEGER i      

C      Integrate[x^2*Exp[-x^2], {x, -\[Infinity], \[Infinity]}]

      anal = SQRT(ATAN(1.0)) ! 0.88622692545275801365
      vsota = 0.0D0
      DO i = 1,GHquad%n 
            vsota = vsota + GHquad%xi(i)*GHquad%xi(i)*GHquad%wi(i)
      END DO

      print *,"abs error = ",ABS(vsota-anal)

      END SUBROUTINE
C
C -----------------------------------------------------------------------------------------
C
      SUBROUTINE GHinit(n)

      INTEGER n

      GHquad%n = n 
      ALLOCATE (GHquad%xi(GHquad%n))
      ALLOCATE (GHquad%wi(GHquad%n))


      select case (GHquad%n)
   
      case (2) 
            GHquad%xi(1)=-0.7071067811865475244008D0
            GHquad%xi(2)= 0.7071067811865475244008D0
            GHquad%wi(1)= 0.8862269254527580136491D0
            GHquad%wi(2)= 0.8862269254527580136491D0
                        
      case (3)
            GHquad%xi(1)=-1.224744871391589049099D0
            GHquad%xi(2)=0.0D0
            GHquad%xi(3)=1.224744871391589049099D0
            GHquad%wi(1)=0.295408975150919337883D0
            GHquad%wi(2)=1.181635900603677351532D0
            GHquad%wi(3)=0.295408975150919337883D0

      case (4)
            GHquad%xi(1)= 1.650680123885784555883D0
            GHquad%xi(2)=-0.5246476232752903178841D0
            GHquad%xi(3)= 0.5246476232752903178841D0
            GHquad%xi(4)= 1.650680123885784555883D0
            GHquad%wi(1)= 0.08131283544724517714304D0
            GHquad%wi(2)= 0.8049140900055128365061D0
            GHquad%wi(3)= 0.8049140900055128365061D0
            GHquad%wi(4)= 0.08131283544724517714304D0         	

      case (5)
            GHquad%xi(1)=-2.020182870456085632929D0  
            GHquad%xi(2)=-0.9585724646138185071128D0  
            GHquad%xi(3)=0.0D0
            GHquad%xi(4)=0.9585724646138185071128D0  
            GHquad%xi(5)=2.020182870456085632929D0  
            GHquad%wi(1)=0.01995324205904591320774D0  
            GHquad%wi(2)=0.3936193231522411598285D0  
            GHquad%wi(3)=0.9453087204829418812257D0  
            GHquad%wi(4)=0.3936193231522411598285D0  
            GHquad%wi(5)=0.01995324205904591320774D0               
           
      case (7)
            GHquad%xi(1)=-2.651961356835233492447D0
            GHquad%xi(2)=-1.673551628767471445032D0
            GHquad%xi(3)=-0.8162878828589646630387D0
            GHquad%xi(4)=0.0D0
            GHquad%xi(5)=0.8162878828589646630387D0
            GHquad%xi(6)=1.673551628767471445032D0
            GHquad%xi(7)=2.651961356835233492447D0

            GHquad%wi(1)=9.71781245099519154149D-4
            GHquad%wi(2)=0.05451558281912703059218D0
            GHquad%wi(3)=0.4256072526101278005203D0
            GHquad%wi(4)=0.810264617556807326765D0
            GHquad%wi(5)=0.4256072526101278005203D0
            GHquad%wi(6)=0.05451558281912703059218D0
            GHquad%wi(7)=9.71781245099519154149D-4

      case (9)
            GHquad%xi(1)=-3.19099320178152760723D0
            GHquad%xi(2)=-2.266580584531843111802D0
            GHquad%xi(3)=-1.468553289216667931667D0
            GHquad%xi(4)=-0.7235510187528375733226D0
            GHquad%xi(5)=0.0D0
            GHquad%xi(6)=0.7235510187528375733226D0
            GHquad%xi(7)=1.468553289216667931667D0
            GHquad%xi(8)=2.266580584531843111802D0
            GHquad%xi(9)=3.19099320178152760723D0

            GHquad%wi(1)=3.960697726326438190459D-5
            GHquad%wi(2)=0.00494362427553694721722D0
            GHquad%wi(3)=0.088474527394376573288D0
            GHquad%wi(4)=0.4326515590025557501998D0
            GHquad%wi(5)=0.7202352156060509571243D0
            GHquad%wi(6)=0.4326515590025557501998D0
            GHquad%wi(7)=0.088474527394376573288D0
            GHquad%wi(8)=0.00494362427553694721722D0        
            GHquad%wi(9)=3.960697726326438190459D-5

      case (11)

            GHquad%xi(1)=-3.668470846559582518458D0
            GHquad%xi(2)=-2.783290099781651770837D0
            GHquad%xi(3)=-2.025948015825755335166D0
            GHquad%xi(4)=-1.32655708449493285595D0
            GHquad%xi(5)=-0.6568095668820997650246D0
            GHquad%xi(6)=0.0D0
            GHquad%xi(7)=0.6568095668820997650246D0  
            GHquad%xi(8)=1.32655708449493285595D0          
            GHquad%xi(9)=2.025948015825755335166D0
            GHquad%xi(10)=2.783290099781651770837D0
            GHquad%xi(11)=3.668470846559582518458D0

            GHquad%wi(1)=1.439560393714258220331D-6
            GHquad%wi(2)=3.468194663233455106434D-4
            GHquad%wi(3)=0.01191139544491153245039D0
            GHquad%wi(4)=0.1172278751677085033818D0
            GHquad%wi(5)=0.4293597523561250284461D0
            GHquad%wi(6)=0.6547592869145917792039D0
            GHquad%wi(7)=0.4293597523561250284461D0  
            GHquad%wi(8)=0.1172278751677085033818D0          
            GHquad%wi(9)=0.01191139544491153245039D0
            GHquad%wi(10)=3.468194663233455106434D-4
            GHquad%wi(11)=1.439560393714258220331D-6

      case default
            print *, "Invalid number of GH points!" 

      end select

      END SUBROUTINE

      END MODULE
