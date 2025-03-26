MODULE reactions

CONTAINS

    !The inputs and outputs of all reactions must be the same. In particular, in each reaction we will have:

    ! Inputs
    ! ======

    !   nh :: Number of altitude levels
    !   p(nh) :: Pressure (Pa)
    !   t(nh) :: Temperature (K)
    !   dens(nh) :: Number density of any ambient species(m-3)


    ! Outputs
    ! =======

    !   rrates(nh) :: Reaction rates in each altitude level (s-1)
    !
    !   rtype :: Reaction type
    !               - rtype = 1
    !                   photodissociations
    !                   or reactions a + c -> b + c
    !                   or reactions a + ice -> b + c
    !               - rtype = 2
    !                   reactions a + a -> b + c
    !               - rtype = 3
    !                   reactions a + b -> c + d
    !
    !   ns :: Number of source or parent species (1 or 2)
    !   sID(2) :: gas ID of the parent species
    !   sISO(2) :: isotope ID of the parent species
    !   npr :: Number of products (1 or 2)
    !   pID(2) :: gas ID of the products
    !   pISO(2) :: isotope ID of the products
    !   pf(2) :: number of molecules of each product species
    !   ref :: Reference for the reaction rate


    !**********************************************************************************************
    !O + O2 + CO2 -> O3 + CO2
    subroutine reaction0001(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + O2 + CO2 -> O3 + CO2
        !==========================================================================================

        !co2/n2 efficiency as a third body = 2.075

        !Inputs
        integer, intent(in) :: nh                           
        real, intent(in) :: t(nh)    
        double precision, intent(in) :: p(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 2.075*6.0e-34*(t(:)/300.)**(-2.4)*dens(:)
        rtype = 3
        
        ns = 2
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 7
        sISO(2) = 0
        sf(2) = 1.0

        npr= 1
        pID(1) = 3
        pISO(1) = 0
        pf(1) = 1.0

        ref = 'sehested et al., j. geophys. res., 100, 1995'


    end subroutine reaction0001

    !**********************************************************************************************
    !O + O + CO2 -> O2 + CO2
    subroutine reaction0002(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + O + CO2 -> O2 + CO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh  
        real, intent(in) :: t(nh)                               
        double precision, intent(in) :: p(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 2.5*9.46e-34*exp(485./t(:))*dens(:) ! nist expression

        rtype = 2

        ns = 1
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 2.0

        npr = 1
        pID(1) = 7
        pISO(1) = 0
        pf(1) = 1.0

        ref = 'NIST kinetics database'


    end subroutine reaction0002

    !**********************************************************************************************
    !O + O3 -> O2 + O2
    subroutine reaction0003(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + O3 -> O2 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh    
        real, intent(in) :: t(nh)                             
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 8.0d-12
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 2060.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 3
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 7
        pISO(1) = 0
        pf(1) = 2.0

        ref = 'JPL 2020'


    end subroutine reaction0003

    !**********************************************************************************************
    !O(1D) + CO2 -> O + CO2
    subroutine reaction0004(nh,p,t,co2,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O(1D) + CO2  -> O + CO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),co2(nh)    

        !Local
        double precision, parameter :: alpha = 7.5d-11
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = -115.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:)) * co2(:)

        rtype = 1

        ns = 1
        sID(1) = 133
        sISO(1) = 0
        sf(1) = 1.0

        npr = 1
        pID(1) = 45
        pISO(1) = 0
        pf(1) = 1.0

        ref = 'JPL 2020'


    end subroutine reaction0004

    !**********************************************************************************************
    !O(1D) + H2O -> OH + OH
    subroutine reaction0005(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O(1D) + H2O  -> OH + OH
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 1.63d-10
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = -60.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 133
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 1
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 13
        pISO(1) = 0
        pf(1) = 2.0

        ref = 'JPL 2020'


    end subroutine reaction0005

    !**********************************************************************************************
    !O(1D) + H2 -> OH + H
    subroutine reaction0006(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O(1D) + H2  -> OH + H
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 1.2d-10
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 0.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 133
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 39
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 13
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 48
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0006

    !**********************************************************************************************
    !O(1D) + O2 -> O + O2
    subroutine reaction0007(nh,p,t,o2,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O(1D) + O2  -> O + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh  
        real, intent(in) :: t(nh)                               
        double precision, intent(in) :: p(nh),o2(nh)    

        !Local
        double precision, parameter :: alpha = 3.3d-11
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = -55.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:)) * o2(:)

        rtype = 1

        ns = 1
        sID(1) = 133
        sISO(1) = 0
        sf(1) = 1.0

        npr = 1
        pID(1) = 45
        pISO(1) = 0
        pf(1) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0007

    !**********************************************************************************************
    !O(1D) + O3  -> O2 + O2
    subroutine reaction0008(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O(1D) + O3  -> O2 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh 
        real, intent(in) :: t(nh)                                
        double precision, intent(in) :: p(nh),dens(nh)  
        
        !Local
        double precision, parameter :: alpha = 2.4d-10
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 0.d0
        double precision, parameter :: br = 0.5d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 133
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 3
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 7
        pISO(1) = 0
        pf(1) = 2.0

        ref = 'JPL 2020'

    end subroutine reaction0008

    !**********************************************************************************************
    !O(1D) + O3  -> O2 + O + O
    subroutine reaction0009(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O(1D) + O3  -> O2 + O + O
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh 
        real, intent(in) :: t(nh)                                
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 2.4d-10
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 0.d0
        double precision, parameter :: br = 0.5d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 133
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 3
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 7
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 45
        pISO(2) = 0
        pf(2) = 2.0

        ref = 'JPL 2020'

    end subroutine reaction0009

    !**********************************************************************************************
    !O + HO2 -> OH + O2
    subroutine reaction0010(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + HO2 -> OH + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)   
        
        !Local
        double precision, parameter :: alpha = 3.0d-11
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = -200.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 44
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 7
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 13
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0010

    !**********************************************************************************************
    !O + OH -> O2 + H
    subroutine reaction0011(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + OH -> O2 + H
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)   
        
        !Local
        double precision, parameter :: alpha = 1.8d-11
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = -180.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 13
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 7
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 48
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0011

    !**********************************************************************************************
    !H + O3 -> OH + O2
    subroutine reaction0012(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       H + O3 -> OH + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh  
        real, intent(in) :: t(nh)                               
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 1.4d-10
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 470.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 48
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 3
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 13
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 7
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0012

    !**********************************************************************************************
    !H + HO2 -> OH + OH
    subroutine reaction0013(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       H + HO2 -> OH + OH
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh    
        real, intent(in) :: t(nh)                             
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 7.2d-11
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 0.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 48
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 44
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 13
        pISO(1) = 0
        pf(1) = 2.0

        ref = 'JPL 2020'

    end subroutine reaction0013

    !**********************************************************************************************
    !H + HO2 -> H2 + O2
    subroutine reaction0014(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       H + HO2 -> H2 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)  
        
        !Local
        double precision, parameter :: alpha = 6.9d-12
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 0.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 48
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 44
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 39
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 7
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0014

    !**********************************************************************************************
    !H + HO2 -> H2O + O
    subroutine reaction0015(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       H + HO2 -> H2O + O
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh  
        real, intent(in) :: t(nh)                               
        double precision, intent(in) :: p(nh),dens(nh)   

        !Local
        double precision, parameter :: alpha = 1.6d-12
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 0.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 48
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 44
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 1
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 45
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0015

    !**********************************************************************************************
    !OH + HO2 -> H2O + O2
    subroutine reaction0016(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       OH + HO2 -> H2O + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 4.8d-11
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = -250.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 13
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 44
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 1
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 7
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0016

    !**********************************************************************************************
    !HO2 + HO2 -> H2O2 + O2
    subroutine reaction0017(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       HO2 + HO2 -> H2O2 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 3.0d-13
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = -460.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 2

        ns = 1
        sID(1) = 44
        sISO(1) = 0
        sf(1) = 2.0

        npr = 2
        pID(1) = 25
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 7
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0017

    !**********************************************************************************************
    !OH + H2O2 -> H2O + HO2
    subroutine reaction0018(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       OH + H2O2 -> H2O + HO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh      
        real, intent(in) :: t(nh)                           
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 1.8d-12
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 0.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 13
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 25
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 1
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 44
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0018

    !**********************************************************************************************
    !OH + H2 -> H2O + H
    subroutine reaction0019(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       OH + H2 -> H2O + H
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh     
        real, intent(in) :: t(nh)                            
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 2.8d-12
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 1800.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 13
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 39
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 1
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 48
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0019

    !**********************************************************************************************
    !H + O2 + CO2 -> HO2 + CO2
    subroutine reaction0020(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       H + O2 + CO2 -> HO2 + CO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        integer :: ih
        double precision :: k0x,kinfx,kf
        double precision, parameter :: k0 = 5.3d-32
        double precision, parameter :: n = 1.8
        double precision, parameter :: kinf = 9.5d-11
        double precision, parameter :: m = -0.4

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        do ih=1,nh

            !k0x = 2.4*4.4e-32*(t(ilev)/300.)**(-1.3) ! FL li et al 2017
            k0x = 2.4*k0*(298./t(ih))**(n)
            kinfx = kinf*(298./t(ih))**(m)
            
            kf = (kinfx*k0x*dens(ih)/(kinfx + k0x*dens(ih)))  &
                *0.6**(1. + (log10(k0x*dens(ih)/kinfx))**2.)**(-1.0)
            
            rrates(ih) = kf
            
        end do

        rtype = 3

        ns = 2
        sID(1) = 48
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 7
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 44
        pISO(1) = 0
        pf(1) = 1.0

        ref = 'JPL 2020'


    end subroutine reaction0020

    !**********************************************************************************************
    !O + H2O2 -> OH + HO2
    subroutine reaction0021(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + H2O2 -> OH + HO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)  
        
        !Local
        double precision, parameter :: alpha = 1.4d-12
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 2000.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 25
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 13
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 44
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0021

    !**********************************************************************************************
    !OH + OH -> H2O + O
    subroutine reaction0022(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       OH + OH -> H2O + O
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 1.8d-12
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 0.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 2

        ns = 1
        sID(1) = 13
        sISO(1) = 0
        sf(1) = 2.0

        npr = 2
        pID(1) = 1
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 45
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0022

    !**********************************************************************************************
    !OH + O3 -> HO2 + O2
    subroutine reaction0023(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       OH + O3 -> HO2 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh  
        real, intent(in) :: t(nh)                               
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 1.7d-12
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 940.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 13
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 3
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 44
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 7
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0023

    !**********************************************************************************************
    !HO2 + O3 -> OH + O2 + O2
    subroutine reaction0024(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       HO2 + O3 -> OH + O2 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh     
        real, intent(in) :: t(nh)                            
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 1.0d-14
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 490.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 44
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 3
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 13
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 7
        pISO(2) = 0
        pf(2) = 2.0

        ref = 'JPL 2020'

    end subroutine reaction0024

    !**********************************************************************************************
    !HO2 + HO2 + CO2 -> H2O2 + O2 + CO2
    subroutine reaction0025(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       HO2 + HO2 + CO2 -> H2O2 + O2 + CO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 2.1d-33
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = -920.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 2.5 * alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:)) * dens(:)

        rtype = 2

        ns = 1
        sID(1) = 44
        sISO(1) = 0
        sf(1) = 2.0

        npr = 2
        pID(1) = 25
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 7
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0025

    !**********************************************************************************************
    !OH + OH + CO2 -> H2O2 + CO2
    subroutine reaction0026(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       OH + OH + CO2 -> H2O2 + CO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        integer :: ih
        double precision :: k0x,kinfx,kf
        double precision, parameter :: k0 = 6.9d-31
        double precision, parameter :: n = 1.0
        double precision, parameter :: kinf = 2.6d-11
        double precision, parameter :: m = 0.0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

         do ih=1,nh

            k0x = 2.5*k0*(298./t(ih))**(n)
            kinfx = kinf*(298./t(ih))**(m)
            
            kf = (kinfx*k0x*dens(ih)/(kinfx + k0x*dens(ih)))  &
                *0.6**(1. + (log10(k0x*dens(ih)/kinfx))**2.)**(-1.0)
            
            rrates(ih) = kf
            
        end do

         rtype = 2

         ns = 1
         sID(1) = 13
         sISO(1) = 0
         sf(1) = 2.0
 
         npr = 1
         pID(1) = 25
         pISO(1) = 0
         pf(1) = 1.0

        ref = 'JPL 2020'
        
    end subroutine reaction0026

    !**********************************************************************************************
    !H + H + CO2 -> H2 + CO2
    subroutine reaction0027(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       H + H + CO2 -> H2 + CO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 2.5*1.8e-30*(t(:)**(-1.0))*dens(:)

        rtype = 2

        ns = 1
        sID(1) = 48
        sISO(1) = 0
        sf(1) = 2.0

        npr = 1
        pID(1) = 39
        pISO(1) = 0
        pf(1) = 1.0

        ref = 'Baulch et al., 2005'

    end subroutine reaction0027

    !**********************************************************************************************
    !O + NO2 + M -> NO + O2 + M
    subroutine reaction0028(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + NO2 + M -> NO + NO + M
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        integer :: ih
        double precision :: k0x,kinfx,kf,kint,kca
        double precision, parameter :: k0 = 3.4d-31
        double precision, parameter :: n = 1.6
        double precision, parameter :: kinf = 2.3d-11
        double precision, parameter :: m = 0.2
        double precision, parameter :: A = 5.3d-12
        double precision, parameter :: B = -200.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        do ih=1,nh

            !association (NO3)
            
            k0x = 2.5*k0*(298./t(ih))**(n)
            kinfx = kinf*(298./t(ih))**(m)
            
            kf = (kinfx*k0x*dens(ih)/(kinfx + k0x*dens(ih)))  &
                *0.6**(1. + (log10(k0x*dens(ih)/kinfx))**2.)**(-1.0)
            
            !chemical activation (NO + O2)
            
            kint = A*exp(-B/t(ih))
            kca = kint*(1.d0 - kf/kinf)
            
            !total : chemical activation
            
            rrates(ih) = kca
            
        end do

        !rrates(:) = 5.1e-12*exp(210./t(:)) #JPL-2006

        rtype = 3

        ns = 2
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 10
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 8
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 7
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0028

    !**********************************************************************************************
    !NO + O3 -> NO2 + O2
    subroutine reaction0029(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       NO + O3 -> NO2 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh         
        real, intent(in) :: t(nh)                        
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 3.0d-12
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 1500.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3
    
        ns = 2
        sID(1) = 8
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 3
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 10
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 7
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2006'

    end subroutine reaction0029

    !**********************************************************************************************
    !NO + HO2 -> NO2 + OH
    subroutine reaction0030(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       NO + HO2 -> NO2 + OH
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 3.44d-12
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = -260.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3
    
        ns = 2
        sID(1) = 8
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 44
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 10
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 13
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2011'

    end subroutine reaction0030

    !**********************************************************************************************
    !N + NO -> N2 + O
    subroutine reaction0031(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       N + NO -> N2 + O
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh  
        real, intent(in) :: t(nh)                               
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 2.1d-11
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = -100.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3
    
        ns = 2
        sID(1) = 47
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 8
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 22
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 45
        pISO(2) = 0
        pf(2) = 1.0


        ref = 'JPL 2020'

    end subroutine reaction0031

    !**********************************************************************************************
    !N + O2 -> NO + O
    subroutine reaction0032(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       N + O2 -> NO + O
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 3.3d-12
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 3150.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 47
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 7
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 8
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 45
        pISO(2) = 0
        pf(2) = 1.0


        ref = 'JPL 2020'

    end subroutine reaction0032

    !**********************************************************************************************
    !NO2 + H -> NO + OH
    subroutine reaction0033(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       NO2 + H -> NO + OH
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh    
        real, intent(in) :: t(nh)                             
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 1.35d-10
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 0.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))
        !rrates(:) = 4.0e-10*exp(-340./t(:)) 

        rtype = 3
    
        ns = 2
        sID(1) = 10
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 48
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 8
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 13
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0033

    !**********************************************************************************************
    !N + O -> NO
    subroutine reaction0034(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       N + O -> NO
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh     
        real, intent(in) :: t(nh)                            
        double precision, intent(in) :: p(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 2.8e-17*(300./t(:))**0.5

        rtype = 3
    
        ns = 2
        sID(1) = 47
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 45
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 8
        pISO(1) = 0
        pf(1) = 1.0


        ref = 'JPL 2011'

    end subroutine reaction0034

    !**********************************************************************************************
    !N + HO2 -> NO + OH
    subroutine reaction0035(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       N + HO2 -> NO + OH
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh       
        real, intent(in) :: t(nh)                          
        double precision, intent(in) :: p(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 2.19e-11

        rtype = 3
    
        ns = 2
        sID(1) = 47
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 44
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 8
        pISO(1) = 0
        pf(1) = 1.0
        pID(1) = 13
        pISO(1) = 0
        pf(1) = 1.0


        ref = 'brune et al., j. chem. phys., 87, 1983'

    end subroutine reaction0035

    !**********************************************************************************************
    !N + OH -> NO + H
    subroutine reaction0036(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       N + OH -> NO + H
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh    
        real, intent(in) :: t(nh)                             
        double precision, intent(in) :: p(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 3.8e-11*exp(85./t(:))

        rtype = 3
    
        ns = 2
        sID(1) = 47
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 13
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 8
        pISO(1) = 0
        pf(1) = 1.0
        pID(1) = 48
        pISO(1) = 0
        pf(1) = 1.0


        ref = 'atkinson et al., j. phys. chem. ref. data, 18, 881, 1989'

    end subroutine reaction0036

    !**********************************************************************************************
    !N(2D) + O  -> N + O
    subroutine reaction0037(nh,p,t,o,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       N(2D) + O  -> N + O
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh       
        real, intent(in) :: t(nh)                          
        double precision, intent(in) :: p(nh),o(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 3.3e-12*exp(-260./t(:)) * o(:)

        rtype = 1
    
        ns = 1
        sID(1) = 134
        sISO(1) = 0
        sf(1) = 1.0

        npr = 1
        pID(1) = 47
        pISO(1) = 0
        pf(1) = 1.0


        ref = 'herron, j. phys. chem. ref. data, 1999'

    end subroutine reaction0037

    !**********************************************************************************************
    !N(2D) + N2  -> N + N2
    subroutine reaction0038(nh,p,t,n2,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       N(2D) + N2  -> N + N2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh 
        real, intent(in) :: t(nh)                                
        double precision, intent(in) :: p(nh),n2(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 1.7e-14 * n2(:)

        rtype = 1
    
        ns = 1
        sID(1) = 134
        sISO(1) = 0
        sf(1) = 1.0

        npr = 1
        pID(1) = 47
        pISO(1) = 0
        pf(1) = 1.0


        ref = 'herron, j. phys. chem. ref. data, 1999'

    end subroutine reaction0038

    !**********************************************************************************************
    !N(2D) + CO2  -> NO + CO
    subroutine reaction0039(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       N(2D) + CO2  -> NO + CO
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh    
        real, intent(in) :: t(nh)                             
        double precision, intent(in) :: p(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 3.6e-13

        rtype = 3
    
        ns = 2
        sID(1) = 134
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 2
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 8
        pISO(1) = 0
        pf(1) = 1.0
        pID(1) = 5
        pISO(1) = 0
        pf(1) = 1.0


        ref = 'herron, j. phys. chem. ref. data, 1999'

    end subroutine reaction0039

    !**********************************************************************************************
    !OH + CO -> CO2 + H
    subroutine reaction0040(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       OH + CO -> CO2 + H
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        integer :: ih
        double precision :: k0x,kinfx,kf,kint,kca
        double precision, parameter :: k0 = 6.9d-33
        double precision, parameter :: n = 2.1
        double precision, parameter :: kinf = 1.1d-12
        double precision, parameter :: m = -1.3
        double precision, parameter :: A = 1.85d-13
        double precision, parameter :: B = 65.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        do ih=1,nh

            !association (HOCO)
            
            k0x = 2.5*k0*(298./t(ih))**(n)
            kinfx = kinf*(298./t(ih))**(m)
            
            kf = (kinfx*k0x*dens(ih)/(kinfx + k0x*dens(ih)))  &
                *0.6**(1. + (log10(k0x*dens(ih)/kinfx))**2.)**(-1.0)
            
            !chemical activation (CO2 + H)
            
            kint = A*exp(-B/t(ih))
            kca = kint*(1.d0 - kf/kinf)
            
            !total : chemical activation
            
            rrates(ih) = kca
            
        end do

        rtype = 3

        ns = 2
        sID(1) = 13
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 5
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 2
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 48
        pISO(2) = 0
        pf(2) = 1.0


        ref = 'JPL 2020'

    end subroutine reaction0040

    !**********************************************************************************************
    !OH + CO -> HOCO
    subroutine reaction0041(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       OH + CO -> HOCO
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        integer :: ih
        double precision :: k0x,kinfx,kf,kint,kca
        double precision, parameter :: k0 = 6.9d-33
        double precision, parameter :: n = 2.1
        double precision, parameter :: kinf = 1.1d-12
        double precision, parameter :: m = -1.3
        double precision, parameter :: A = 1.85d-13
        double precision, parameter :: B = 65.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        do ih=1,nh

            !association (HOCO)
            
            k0x = 2.5*k0*(298./t(ih))**(n)
            kinfx = kinf*(298./t(ih))**(m)
            
            kf = (kinfx*k0x*dens(ih)/(kinfx + k0x*dens(ih)))  &
                *0.6**(1. + (log10(k0x*dens(ih)/kinfx))**2.)**(-1.0)
            
            !chemical activation (CO2 + H)
            
            kint = A*exp(-B/t(ih))
            kca = kint*(1.d0 - kf/kinf)
            
            !total : chemical activation
            
            rrates(ih) = kf
            
        end do

        rtype = 3

        ns = 2
        sID(1) = 13
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 5
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 80
        pISO(1) = 0
        pf(1) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0041

    !**********************************************************************************************
    !O + CO + M -> CO2 + M
    subroutine reaction0042(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + CO + M -> CO2 + M
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh    
        real, intent(in) :: t(nh)                             
        double precision, intent(in) :: p(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 2.5*6.5e-33*exp(-2184./t(:))*dens(:)

        rtype = 3

        ns = 2
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 5
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 2
        pISO(1) = 0
        pf(1) = 1.0


        ref = 'tsang and hampson, 1986'

    end subroutine reaction0042

    !**********************************************************************************************
    !O(1D) + N2 + CO2 -> N2O + CO2
    subroutine reaction0043(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O(1D) + N2 + M -> N2O + M
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        integer :: ih
        double precision :: k0x,kinfx,kf,kint,kca
        double precision, parameter :: k0 = 2.8d-36
        double precision, parameter :: n = 0.9
        double precision, parameter :: kinf = 0.d0
        double precision, parameter :: m = 0.d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 2.5*k0*(t(:)/300.)**(-n)*dens(:)

        rtype = 3

        ns = 2
        sID(1) = 133
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 22
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 4
        pISO(1) = 0
        pf(1) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0043

    !**********************************************************************************************
    !O + NO + CO2 -> NO2 + CO2
    subroutine reaction0044(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + NO + CO2 -> NO2 + CO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        integer :: ih
        double precision :: k0x,kinfx,kf
        double precision, parameter :: k0 = 9.1d-32
        double precision, parameter :: n = 1.5
        double precision, parameter :: kinf = 3.0d-11
        double precision, parameter :: m = 0.0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        do ih=1,nh

            k0x = 2.4*k0*(298./t(ih))**(n)
            kinfx = kinf*(298./t(ih))**(m)
            
            kf = (kinfx*k0x*dens(ih)/(kinfx + k0x*dens(ih)))  &
                *0.6**(1. + (log10(k0x*dens(ih)/kinfx))**2.)**(-1.0)
            
            rrates(ih) = kf
            
        end do

        rtype = 3

        ns = 2
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 8
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 10
        pISO(1) = 0
        pf(1) = 1.0

        ref = 'JPL 2020'


    end subroutine reaction0044

    !**********************************************************************************************
    !O(1D) + N2 -> O + N2
    subroutine reaction0045(nh,p,t,n2,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O(1D) + N2 -> O(1D) + N2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),n2(nh)    

        !Local
        double precision, parameter :: alpha = 2.5d-11
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = -110.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:)) * n2(:)

        rtype = 1

        ns = 1
        sID(1) = 133
        sISO(1) = 0
        sf(1) = 1.0

        npr = 1
        pID(1) = 45
        pISO(1) = 0
        pf(1) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0045

    !**********************************************************************************************
    !O(1D) + N2O -> N2 + O2
    subroutine reaction0046(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O(1D) + N2O -> N2 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 1.19d-10
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = -20.d0
        double precision, parameter :: br = 0.39d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 133
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 4
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 22
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 7
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0046

    !**********************************************************************************************
    !O(1D) + N2O -> NO + NO
    subroutine reaction0047(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O(1D) + N2O -> NO + NO
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 1.19d-10
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = -20.d0
        double precision, parameter :: br = 0.61d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 133
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 4
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 8
        pISO(1) = 0
        pf(1) = 2.0

        ref = 'JPL 2020'

    end subroutine reaction0047

    !**********************************************************************************************
    !O + NO2 + M -> NO + O2 + M
    subroutine reaction0048(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + NO2 + M -> NO + NO + M
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        integer :: ih
        double precision :: k0x,kinfx,kf,kint,kca
        double precision, parameter :: k0 = 3.4d-31
        double precision, parameter :: n = 1.6
        double precision, parameter :: kinf = 2.3d-11
        double precision, parameter :: m = 0.2
        double precision, parameter :: A = 5.3d-12
        double precision, parameter :: B = -200.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        do ih=1,nh

            !association (NO3)
            
            k0x = 2.5*k0*(298./t(ih))**(n)
            kinfx = kinf*(298./t(ih))**(m)
            
            kf = (kinfx*k0x*dens(ih)/(kinfx + k0x*dens(ih)))  &
                *0.6**(1. + (log10(k0x*dens(ih)/kinfx))**2.)**(-1.0)
            
            !chemical activation (NO + O2)
            
            kint = A*exp(-B/t(ih))
            kca = kint*(1.d0 - kf/kinf)
            
            !total : chemical activation
            
            rrates(ih) = kca
            
        end do

        rtype = 3

        ns = 2
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 10
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 8
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 7
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0048

    !**********************************************************************************************
    !O + NO2 + M -> NO3 + M
    subroutine reaction0049(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + NO2 + M -> NO3 + M
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        integer :: ih
        double precision :: k0x,kinfx,kf,kint,kca
        double precision, parameter :: k0 = 3.4d-31
        double precision, parameter :: n = 1.6
        double precision, parameter :: kinf = 2.3d-11
        double precision, parameter :: m = 0.2
        double precision, parameter :: A = 5.3d-12
        double precision, parameter :: B = -200.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        do ih=1,nh

            !association (NO3)
            
            k0x = 2.5*k0*(298./t(ih))**(n)
            kinfx = kinf*(298./t(ih))**(m)
            
            kf = (kinfx*k0x*dens(ih)/(kinfx + k0x*dens(ih)))  &
                *0.6**(1. + (log10(k0x*dens(ih)/kinfx))**2.)**(-1.0)
            
            !chemical activation (NO + O2)
            
            kint = A*exp(-B/t(ih))
            kca = kint*(1.d0 - kf/kinf)
            
            !total : association
            
            rrates(ih) = kf
            
        end do

        rtype = 3

        ns = 2
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 10
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 91
        pISO(1) = 0
        pf(1) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0049

    !**********************************************************************************************
    !O + NO3 -> O2 + NO2
    subroutine reaction0050(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + NO3 -> O2 + NO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 1.3d-11
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 0.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 91
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 8
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 10
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0050

    !**********************************************************************************************
    !N + NO2 -> N2O + O
    subroutine reaction0051(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       N + NO2 -> N2O + O
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 5.8d-12
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = -220.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 47
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 10
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 4
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 45
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0051

    !**********************************************************************************************
    !NO + NO3 -> NO2 + NO2
    subroutine reaction0052(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       NO + NO3 -> NO2 + NO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 1.7d-11
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = -125.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 8
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 91
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 10
        pISO(1) = 0
        pf(1) = 2.0

        ref = 'JPL 2020'

    end subroutine reaction0052

    !**********************************************************************************************
    !NO2 + O3 -> NO3 + O2
    subroutine reaction0053(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       NO2 + O3 -> NO3 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 1.2d-13
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 2450.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 10
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 3
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 91
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 7
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0053

    !**********************************************************************************************
    !NO3 + NO3 -> 2NO2 + O2
    subroutine reaction0054(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       NO3 + NO3 -> 2NO2 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 8.5d-13
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 2450.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 2

        ns = 1
        sID(1) = 91
        sISO(1) = 0
        sf(1) = 2.0

        npr = 2
        pID(1) = 10
        pISO(1) = 0
        pf(1) = 2.0
        pID(2) = 7
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0054

    !**********************************************************************************************
    !O2 + HOCO -> HO2 + CO2
    subroutine reaction0055(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O2 + HOCO -> HO2 + CO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        double precision, parameter :: alpha = 2.0d-12
        double precision, parameter :: beta = 0.d0
        double precision, parameter :: gamma = 0.d0
        double precision, parameter :: br = 1.0d0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * br * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 2
        sID(1) = 7
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 80
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 44
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 2
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2020'

    end subroutine reaction0055



    !**********************************************************************************************

    !**********************************************************************************************
    !**********************************************************************************************
    !**********************************************************************************************
    !                                      13C isotopes
    !**********************************************************************************************
    !**********************************************************************************************
    !**********************************************************************************************

    !N(2D) + (13C)O2  -> NO + CO
    subroutine reaction0039_13c(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       N(2D) + (13C)O2  -> NO + CO
        !==========================================================================================

        !we assume it is the same rate as reaction0039

        !Inputs
        integer, intent(in) :: nh       
        real, intent(in) :: t(nh)                          
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        integer :: rtype1,ns1,npr1,sID1(2),sISO1(2),pID1(2),pISO1(2)
        double precision :: rrates1(nh)
        real :: sf1(2),pf1(2)
        character (len = 100) :: ref1

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        !Calling reaction r0041 to get reaction rate
        call reaction0039(nh,p,t,dens,rrates1,rtype1,ns1,sID1,sISO1,sf1,npr1,pID1,pISO1,pf1,ref1)

        rrates(:) = rrates1(:)*1.0

        rtype = 3
        
        ns = 2
        sID(1) = 134
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 2
        sISO(2) = 2
        sf(2) = 1.0

        npr = 2
        pID(1) = 8
        pISO(1) = 0
        pf(1) = 1.0
        pID(1) = 5
        pISO(1) = 2
        pf(1) = 1.0

        ref = 'N/A'

    end subroutine reaction0039_13c

    !**********************************************************************************************
    !OH + (13C)O -> (13C)O2 + H
    subroutine reaction0040_13c(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       OH + (13C)O -> (13C)O2 + H
        !==========================================================================================

        !Assumed to be the same as e001, but with a 
        !fractionation factor from Stevens et al. (1980)

        !A polynomial function is fit to capture the pressure-
        !dependence of the fractionation

        !k13 / k12 = 1.00638 - 1.693e-5*press(hPa) + 4.6968e-9 * press(hPa)**2.
        !k13 / k12 = 1.00638 - 1.693e-7*press(Pa) + 4.6968e-13 * press(Pa)**2.

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        integer :: rtype1,ns1,npr1,sID1(2),sISO1(2),pID1(2),pISO1(2)
        double precision :: rrates1(nh)
        real :: sf1(2),pf1(2)
        character (len = 100) :: ref1

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        !Calling reaction r0040 to get reaction rate for oh + co -> co2 + h
        call reaction0040(nh,p,t,dens,rrates1,rtype1,ns1,sID1,sISO1,sf1,npr1,pID1,pISO1,pf1,ref1)

        rrates(:) = rrates1(:)*(1.00638 - 1.693e-7*p(:) + 4.6968e-13*p(:)**2.)

        rtype = 3
    
        ns = 2
        sID(1) = 13
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 5
        sISO(2) = 2
        sf(2) = 1.0

        npr = 2
        pID(1) = 2
        pISO(1) = 2
        pf(1) = 1.0
        pID(2) = 48
        pISO(2) = 0
        pf(2) = 1.0


        ref = 'Stevens et al. (1980)'

    end subroutine reaction0040_13c

    !**********************************************************************************************
    !OH + (13C)O -> HO(13C)O
    subroutine reaction0041_13c(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       OH + (13C)O -> (13C)O2 + H
        !==========================================================================================

        !Assumed to be the same as e001, but with a 
        !fractionation factor from Stevens et al. (1980)

        !A polynomial function is fit to capture the pressure-
        !dependence of the fractionation

        !k13 / k12 = 1.00638 - 1.693e-5*press(hPa) + 4.6968e-9 * press(hPa)**2.
        !k13 / k12 = 1.00638 - 1.693e-7*press(Pa) + 4.6968e-13 * press(Pa)**2.

        !Inputs
        integer, intent(in) :: nh   
        real, intent(in) :: t(nh)                              
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        integer :: rtype1,ns1,npr1,sID1(2),sISO1(2),pID1(2),pISO1(2)
        double precision :: rrates1(nh)
        real :: sf1(2),pf1(2)
        character (len = 100) :: ref1

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        !Calling reaction r0041 to get reaction rate for oh + co -> hoco
        call reaction0041(nh,p,t,dens,rrates1,rtype1,ns1,sID1,sISO1,sf1,npr1,pID1,pISO1,pf1,ref1)

        rrates(:) = rrates1(:)*(1.00638 - 1.693e-7*p(:) + 4.6968e-13*p(:)**2.)

        rtype = 3
    
        ns = 2
        sID(1) = 13
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 5
        sISO(2) = 2
        sf(2) = 1.0

        npr = 2
        pID(1) = 2
        pISO(1) = 2
        pf(1) = 1.0
        pID(2) = 48
        pISO(2) = 0
        pf(2) = 1.0


        ref = 'Stevens et al. (1980)'

    end subroutine reaction0041_13c


    !**********************************************************************************************
    !O + (13C)O + M -> (13C)O2 + M
    subroutine reaction0042_13c(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + (13C)O + M -> (13C)O2 + M
        !==========================================================================================

        !we assume it is the same rate as reaction0041

        !Inputs
        integer, intent(in) :: nh 
        real, intent(in) :: t(nh)                                
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        integer :: rtype1,ns1,npr1,sID1(2),sISO1(2),pID1(2),pISO1(2)
        double precision :: rrates1(nh)
        real :: sf1(2),pf1(2)
        character (len = 100) :: ref1

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        !Calling reaction r0042 to get reaction rate
        call reaction0042(nh,p,t,dens,rrates1,rtype1,ns1,sID1,sISO1,sf1,npr1,pID1,pISO1,pf1,ref1)

        rrates(:) = rrates1(:)*1.0

        rtype = 3
    
        ns = 2
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 5
        sISO(2) = 2
        sf(2) = 1.0

        npr = 1
        pID(1) = 2
        pISO(1) = 2
        pf(1) = 1.0

        ref = 'N/A'

    end subroutine reaction0042_13c




    !**********************************************************************************************
    !**********************************************************************************************
    !**********************************************************************************************
    !                                      15N isotopes
    !**********************************************************************************************
    !**********************************************************************************************
    !**********************************************************************************************

    !**********************************************************************************************
    !(15N)O + O3 --> (15N)O2 + O2
    subroutine reaction0029_15n(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       (15N)O + O3 --> (15N)O2 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh           
        real, intent(in) :: t(nh)                    
        double precision, intent(in) :: p(nh),dens(nh)    

        !Local
        integer :: rtype1,ns1,npr1,sID1(2),sISO1(2),pID1(2),pISO1(2)
        double precision :: rrates1(nh)
        real :: sf1(2),pf1(2)
        character (len = 100) :: ref1

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        real, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        !Calling reaction r0029 to get reaction rate
        call reaction0029(nh,p,t,dens,rrates1,rtype1,ns1,sID1,sISO1,sf1,npr1,pID1,pISO1,pf1,ref1)

        rrates(:) = rrates1(:)*1.0

        rtype = 3

        ns = 3
        sID(1) = 8
        sISO(1) = 2
        sf(1) = 1.0
        sID(2) = 3
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 10
        pISO(1) = 2
        pf(1) = 1.0
        pID(2) = 7
        pISO(2) = 0
        pf(2) = 1.0


        ref = 'Walters and Michalski (2016)'

    end subroutine reaction0029_15n


    
END MODULE reactions