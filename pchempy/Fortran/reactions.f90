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

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
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
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
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
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 8.0e-12*exp(-2060./t(:))

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

        ref = 'JPL 2003'


    end subroutine reaction0003

    !**********************************************************************************************
    !O(1D) + CO2 -> O + CO2
    subroutine reaction0004(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O(1D) + CO2  -> O + CO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 7.5e-11*exp(115./t(:)) * dens(:)

        rtype = 1

        ns = 1
        sID(1) = 133
        sISO(1) = 0
        sf(1) = 1.0

        npr = 1
        pID(1) = 45
        pISO(1) = 0
        pf(1) = 1.0

        ref = 'JPL 2006'


    end subroutine reaction0004

    !**********************************************************************************************
    !O(1D) + H2O -> OH + OH
    subroutine reaction0005(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O(1D) + H2O  -> OH + OH
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 1.63e-10*exp(60./t(:))

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

        ref = 'JPL 2006'


    end subroutine reaction0005

    !**********************************************************************************************
    !O(1D) + H2 -> OH + H
    subroutine reaction0006(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O(1D) + H2  -> OH + H
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 1.2e-10

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

        ref = 'JPL 2011'

    end subroutine reaction0006

    !**********************************************************************************************
    !O(1D) + O2 -> O + O2
    subroutine reaction0007(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O(1D) + O2  -> O + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 3.3e-11*exp(55./t(:)) * dens(:)

        rtype = 1

        ns = 1
        sID(1) = 133
        sISO(1) = 0
        sf(1) = 1.0

        npr = 1
        pID(1) = 45
        pISO(1) = 0
        pf(1) = 1.0

        ref = 'JPL 2006'

    end subroutine reaction0007

    !**********************************************************************************************
    !O(1D) + O3  -> O2 + O2
    subroutine reaction0008(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O(1D) + O3  -> O2 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 1.2e-10

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

        ref = 'JPL 2003'

    end subroutine reaction0008

    !**********************************************************************************************
    !O(1D) + O3  -> O2 + O + O
    subroutine reaction0009(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O(1D) + O3  -> O2 + O + O
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 1.2e-10

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

        ref = 'JPL 2003'

    end subroutine reaction0009

    !**********************************************************************************************
    !O + HO2 -> OH + O2
    subroutine reaction0010(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + HO2 -> OH + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 3.0e-11*exp(200./t(:))

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

        ref = 'JPL 2003'

    end subroutine reaction0010

    !**********************************************************************************************
    !O + OH -> O2 + H
    subroutine reaction0011(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + OH -> O2 + H
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 1.8e-11*exp(180./t(:))

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

        ref = 'JPL 2011'

    end subroutine reaction0011

    !**********************************************************************************************
    !H + O3 -> OH + O2
    subroutine reaction0012(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       H + O3 -> OH + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 1.4e-10*exp(-470./t(:))

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

        ref = 'JPL 2003'

    end subroutine reaction0012

    !**********************************************************************************************
    !H + HO2 -> OH + OH
    subroutine reaction0013(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       H + HO2 -> OH + OH
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 7.2e-11

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

        ref = 'JPL 2006'

    end subroutine reaction0013

    !**********************************************************************************************
    !H + HO2 -> H2 + O2
    subroutine reaction0014(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       H + HO2 -> H2 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 6.9e-12

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

        ref = 'JPL 2006'

    end subroutine reaction0014

    !**********************************************************************************************
    !H + HO2 -> H2O + O
    subroutine reaction0015(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       H + HO2 -> H2O + O
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 1.6e-12

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

        ref = 'JPL 2006'

    end subroutine reaction0015

    !**********************************************************************************************
    !OH + HO2 -> H2O + O2
    subroutine reaction0016(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       OH + HO2 -> H2O + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 4.8e-11*exp(250./t(:))

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

        ref = 'JPL 2003'

    end subroutine reaction0016

    !**********************************************************************************************
    !HO2 + HO2 -> H2O2 + O2
    subroutine reaction0017(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       HO2 + HO2 -> H2O2 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 3.0e-13*exp(460./t(:))

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

        ref = 'JPL 2015'

    end subroutine reaction0017

    !**********************************************************************************************
    !OH + H2O2 -> H2O + HO2
    subroutine reaction0018(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       OH + H2O2 -> H2O + HO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 1.8e-12

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

        ref = 'JPL 2006'

    end subroutine reaction0018

    !**********************************************************************************************
    !OH + H2 -> H2O + H
    subroutine reaction0019(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       OH + H2 -> H2O + H
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = 2.8e-12*exp(-1800./t(:))

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

        ref = 'JPL 2006'

    end subroutine reaction0019

    !**********************************************************************************************
    !H + O2 + CO2 -> HO2 + CO2
    subroutine reaction0020(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       H + O2 + CO2 -> HO2 + CO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Local
        integer :: ih
        double precision :: ak0,ak1,xpo

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        do ih = 1,nh
            !ak0 = 3.1*2.4*4.4e-32*(t(ilev)/300.)**(-1.3) ! FL li et al 2017
            ak0 = 2.4*4.4e-32*(t(ih)/300.)**(-1.3)
            ak1 = 7.5e-11*(t(ih)/300.)**(0.2)
            
            rate = (ak0*dens(ih))/(1. + ak0*dens(ih)/ak1)
            xpo = 1./(1. + dlog10((ak0*dens(ih))/ak1)**2)
            rrates(ih) = rate*0.6**xpo
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

        ref = 'JPL 2011'

    end subroutine reaction0020

    !**********************************************************************************************
    !O + H2O2 -> OH + HO2
    subroutine reaction0021(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + H2O2 -> OH + HO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 1.4e-12*exp(-2000./t(:))

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

        ref = 'JPL 2003'

    end subroutine reaction0021

    !**********************************************************************************************
    !OH + OH -> H2O + O
    subroutine reaction0022(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       OH + OH -> H2O + O
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 1.8e-12

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

        ref = 'JPL 2006'

    end subroutine reaction0022

    !**********************************************************************************************
    !OH + O3 -> HO2 + O2
    subroutine reaction0023(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       OH + O3 -> HO2 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 1.7e-12*exp(-940./t(:))

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

        ref = 'JPL 2003'

    end subroutine reaction0023

    !**********************************************************************************************
    !HO2 + O3 -> OH + O2 + O2
    subroutine reaction0024(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       HO2 + O3 -> OH + O2 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 1.0e-14*exp(-490./t(:))

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

        ref = 'JPL 2003'

    end subroutine reaction0024

    !**********************************************************************************************
    !HO2 + HO2 + CO2 -> H2O2 + O2 + CO2
    subroutine reaction0025(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       HO2 + HO2 + CO2 -> H2O2 + O2 + CO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 2.5*2.1e-33*exp(920./t(:))*dens(:)

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

        ref = 'JPL 2011'

    end subroutine reaction0025

    !**********************************************************************************************
    !OH + OH + CO2 -> H2O2 + CO2
    subroutine reaction0026(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       OH + OH + CO2 -> H2O2 + CO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Local
        double precision :: ak0,ak1,rate,xpo
        integer :: ih

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        do ih = 1,nh
            ak0 = 2.5*6.9e-31*(t(ih)/300.)**(-1.0)
            ak1 = 2.6e-11*(t(ih)/300.)**(0.0)
   
            rate = (ak0*dens(ih))/(1. + ak0*dens(ih)/ak1)
            xpo = 1./(1. + dlog10((ak0*dens(ih))/ak1)**2)
            rrates(ih) = rate*0.6**xpo
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

        ref = 'JPL 2003'

    end subroutine reaction0026

    !**********************************************************************************************
    !H + H + CO2 -> H2 + CO2
    subroutine reaction0027(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       H + H + CO2 -> H2 + CO2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
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
    !NO2 + O -> NO + O2
    subroutine reaction0028(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       NO2 + O -> NO + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 5.1e-12*exp(210./t(:))

        rtype = 3
    
        ns = 2
        sID(1) = 10
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 45
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 8
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 7
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'JPL 2006'

    end subroutine reaction0028

    !**********************************************************************************************
    !NO + O3 -> NO2 + O2
    subroutine reaction0029(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       NO + O3 -> NO2 + O2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 3.0e-12*exp(-1500./t(:))

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
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 3.3e-12*exp(270./t(:))

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
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 2.1e-11*exp(100./t(:))

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


        ref = 'JPL 2011'

    end subroutine reaction0031

    !**********************************************************************************************
    !N + O2 -> NO + O
    subroutine reaction0032(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       N + O2 -> NO + O
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 1.5e-11*exp(-3600./t(:))

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


        ref = 'JPL 2011'

    end subroutine reaction0032

    !**********************************************************************************************
    !NO2 + H -> NO + OH
    subroutine reaction0033(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       NO2 + H -> NO + OH
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 4.0e-10*exp(-340./t(:))

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


        ref = 'JPL 2011'

    end subroutine reaction0033

    !**********************************************************************************************
    !N + O -> NO
    subroutine reaction0034(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       N + O -> NO
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
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
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
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
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
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
        double precision, intent(in) :: p(nh),t(nh),o(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
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
        double precision, intent(in) :: p(nh),t(nh),n2(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
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
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
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
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Local
        integer :: ih
        double precision :: k0,kinf,kf,kint,kca

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        do ih=1,nh

            !association
            
            k0 = 2.5*6.9d-33*(298./t(ih))**(2.1)
            kinf = 1.1d-12*(298./t(ih))**(-1.3)
            
            kf = (kinf*k0*dens(ih)/(kinf + k0*dens(ih)))  &
                *0.6**(1. + (log10(k0*dens(ih)/kinf))**2.)**(-1.0)
            
            !chemical activation
            
            kint = 1.85e-13*exp(-65./t(ih))
            
            kca = kint*(1.d0 - kf/kinf)
            
            !total : association + chemical activation
            
            rrates(ih) = kf + kca
            
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


        ref = 'JPL 2019'

    end subroutine reaction0040

    !**********************************************************************************************
    !O + CO + M -> CO2 + M
    subroutine reaction0041(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + CO + M -> CO2 + M
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
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

    end subroutine reaction0041

    !**********************************************************************************************
    !C + H -> CH
    subroutine reaction0042(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       C + H -> CH
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Local
        double precision,  parameter :: alpha = 1.0d-17
        double precision,  parameter :: beta = 0.d0
        double precision,  parameter :: gamma = 0.d0


        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = alpha * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 3
        sID(1) = 46
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 48
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 93
        pISO(1) = 0
        pf(1) = 1.0


        ref = 'UMIST RATE12'

    end subroutine reaction0042

    !**********************************************************************************************
    !C + N -> CN
    subroutine reaction0043(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       C + N -> CN
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Local
        double precision,  parameter :: alpha = 5.72d-19
        double precision,  parameter :: beta = 0.37
        double precision,  parameter :: gamma = 51.

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = alpha * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))
    

        rtype = 3

        ns = 3
        sID(1) = 46
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 47
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 96
        pISO(1) = 0
        pf(1) = 1.0


        ref = 'UMIST RATE12'

    end subroutine reaction0043

    !**********************************************************************************************
    !CH + H → H2 + C
    subroutine reaction0044(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        CH + H → H2 + C
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Local
        double precision,  parameter :: alpha = 1.31e-10
        double precision,  parameter :: beta = 0.
        double precision,  parameter :: gamma = 80.

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = alpha * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 3
        sID(1) = 93
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 48
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 39
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 46
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0044 

    !**********************************************************************************************
    !CH + N → CN + H
    subroutine reaction0045(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        CH + N → CN + H
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Local
        double precision,  parameter :: alpha = 1.66d-10
        double precision,  parameter :: beta = -0.09
        double precision,  parameter :: gamma = 0.0

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = alpha * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 3
        sID(1) = 93
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 47
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 96
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 48
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0045

    !**********************************************************************************************
    !CH + O → CO + H
    subroutine reaction0046(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        CH + O → CO + H
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Local
        double precision,  parameter :: alpha = 6.02d-11
        double precision,  parameter :: beta = 0.10
        double precision,  parameter :: gamma = -4.5

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = alpha * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 3
        sID(1) = 93
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 45
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 5
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 48
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0046

    !**********************************************************************************************

    !CH + O → OH + C
    subroutine reaction0047(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        CH + O → OH + C
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Local
        double precision,  parameter :: alpha = 2.52d-11
        double precision,  parameter :: beta = 0.0
        double precision,  parameter :: gamma = 2381.

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        rrates(:) = alpha * ((t(:)/300.d0)**beta) * dexp(-gamma/t(:))

        rtype = 3

        ns = 3
        sID(1) = 93
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 45
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 13
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 46
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0047

    !**********************************************************************************************

    !NH + CN -> HCN + N
    subroutine reaction0048(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        NH + CN -> HCN + N
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 2.94d-12 * (t(:)**(0.5)) * exp(-1000./t(:))

        rtype = 3

        ns = 2
        sID(1) = 113
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 96
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 23
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 47
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0048

    !**********************************************************************************************

    !O + CN → CO + N
    subroutine reaction0049(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        O + CN → CO + N
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 5.00d-11 * (t(:)**(0.0)) * exp(-200./t(:))

        rtype = 3

        ns = 3
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 96
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 5
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 47
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0049

    !**********************************************************************************************

    !O + CN → NO + C
    subroutine reaction0050(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        O + CN → NO + C
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 5.00d-11 * (t(:)**(0.0)) * exp(-200./t(:))

        rtype = 3

        ns = 3
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 96
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 8
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 46
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0050

    !**********************************************************************************************

    !CH + CO2 → HCO + CO
    subroutine reaction0051(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        CH + CO2 → HCO + CO
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 2.94d-13 * (t(:)**(0.5)) * exp(-3000./t(:))

        rtype = 3

        ns = 3
        sID(1) = 93
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 2
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 81
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 5
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0051

    !**********************************************************************************************

    !H + CO2 → CO + OH
    subroutine reaction0052(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        H + CO2 → CO + OH
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 3.38d-10 * (t(:)**(0.0)) * exp(-13163./t(:))

        rtype = 3

        ns = 3
        sID(1) = 48
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 2
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 5
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 13
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0052

    !**********************************************************************************************

    !N + CO2 → NO + CO
    subroutine reaction0053(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        N + CO2 → NO + CO
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 3.2d-13 * (t(:)**(0.0)) * exp(-1710./t(:))

        rtype = 3

        ns = 3
        sID(1) = 47
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 2
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 8
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 5
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0053

    !**********************************************************************************************

    !O + CO2 → O2 + CO
    subroutine reaction0054(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        O + CO2 → O2 + CO
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 2.46d-11 * (t(:)**(0.0)) * exp(-26567./t(:))

        rtype = 3

        ns = 3
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 2
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 7
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 5
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0054

    !**********************************************************************************************

    !H2 + C → CH + H
    subroutine reaction0055(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        H2 + C → CH + H
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 6.64d-10 * (t(:)**(0.0)) * exp(-11700./t(:))

        rtype = 3

        ns = 3
        sID(1) = 39
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 46
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 93
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 48
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0055

    !**********************************************************************************************

    !H2 + O → OH + H
    subroutine reaction0056(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        H2 + O → OH + H
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 3.14d-13 * (t(:)**(2.70)) * exp(-3150./t(:))

        rtype = 3

        ns = 3
        sID(1) = 39
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 45
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 13
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 48
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0056

    !**********************************************************************************************

    !H + H2O → OH + H2
    subroutine reaction0057(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        H + H2O → OH + H2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 1.59d-11 * (t(:)**(1.20)) * exp(-9610./t(:))

        rtype = 3

        ns = 3
        sID(1) = 48
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 1
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 13
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 39
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0057

    !**********************************************************************************************

    !O + H2O → OH + OH
    subroutine reaction0058(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        O + H2O → OH + OH
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 1.85d-11 * (t(:)**(0.95)) * exp(-8571./t(:))

        rtype = 3

        ns = 3
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 1
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 13
        pISO(1) = 0
        pf(1) = 2.0

        ref = 'UMIST RATE12'

    end subroutine reaction0058

    !**********************************************************************************************

    !H + HCN → CN + H2
    subroutine reaction0059(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        H + HCN → CN + H2
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 6.2d-10 * (t(:)**(0.0)) * exp(-12500./t(:))

        rtype = 3

        ns = 3
        sID(1) = 48
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 23
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 96
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 39
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0059

    !**********************************************************************************************

    !O + HCN → CN + OH
    subroutine reaction0060(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        O + HCN → CN + OH
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 6.21d-10 * (t(:)**(0.0)) * exp(-12439./t(:))

        rtype = 3

        ns = 3
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 23
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 96
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 13
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0060

    !**********************************************************************************************

    !O + HCN → CO + NH
    subroutine reaction0061(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !        O + HCN → CO + NH
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        rrates(:) = 7.30d-13 * (t(:)**(1.14)) * exp(-3742./t(:))

        rtype = 3

        ns = 3
        sID(1) = 45
        sISO(1) = 0
        sf(1) = 1.0
        sID(2) = 23
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 5
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 113
        pISO(2) = 0
        pf(2) = 1.0

        ref = 'UMIST RATE12'

    end subroutine reaction0061

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
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Local
        integer :: rtype1,ns1,npr1,sID1(2),sISO1(2),pID1(2),pISO1(2)
        double precision :: sf1(2),pf1(2),rrates1(nh)
        character (len = 100) :: ref1

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
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
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Local
        integer :: rtype1,ns1,npr1,sID1(2),sISO1(2),pID1(2),pISO1(2)
        double precision :: sf1(2),pf1(2),rrates1(nh)
        character (len = 100) :: ref1

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
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
    !O + (13C)O + M -> (13C)O2 + M
    subroutine reaction0041_13c(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       O + (13C)O + M -> (13C)O2 + M
        !==========================================================================================

        !we assume it is the same rate as reaction0041

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Local
        integer :: rtype1,ns1,npr1,sID1(2),sISO1(2),pID1(2),pISO1(2)
        double precision :: sf1(2),pf1(2),rrates1(nh)
        character (len = 100) :: ref1

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref


        !Calling reaction r0041 to get reaction rate
        call reaction0041(nh,p,t,dens,rrates1,rtype1,ns1,sID1,sISO1,sf1,npr1,pID1,pISO1,pf1,ref1)

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

    end subroutine reaction0041_13c

    !**********************************************************************************************
    !(13C) + H -> (13C)H
    subroutine reaction0042_13c(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       (13C) + H -> (13C)H
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Local
        integer :: rtype1,ns1,npr1,sID1(2),sISO1(2),pID1(2),pISO1(2)
        double precision :: sf1(2),pf1(2),rrates1(nh)
        character (len = 100) :: ref1

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        !Calling reaction r0042 to get reaction rate
        call reaction0042(nh,p,t,dens,rrates1,rtype1,ns1,sID1,sISO1,sf1,npr1,pID1,pISO1,pf1,ref1)

        rrates(:) = rrates1(:)*1.0

        rtype = 3

        ns = 3
        sID(1) = 46
        sISO(1) = 2
        sf(1) = 1.0
        sID(2) = 48
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 93
        pISO(1) = 3
        pf(1) = 1.0


        ref = 'N/A'

    end subroutine reaction0042_13c

    !**********************************************************************************************
    !(13C) + N -> (13C)N
    subroutine reaction0043_13c(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       (13C) + N -> (13C)N
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Local
        integer :: rtype1,ns1,npr1,sID1(2),sISO1(2),pID1(2),pISO1(2)
        double precision :: sf1(2),pf1(2),rrates1(nh)
        character (len = 100) :: ref1

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        !Calling reaction r0043 to get reaction rate
        call reaction0043(nh,p,t,dens,rrates1,rtype1,ns1,sID1,sISO1,sf1,npr1,pID1,pISO1,pf1,ref1)

        rrates(:) = rrates1(:)*1.0

        rtype = 3

        ns = 3
        sID(1) = 46
        sISO(1) = 2
        sf(1) = 1.0
        sID(2) = 47
        sISO(2) = 0
        sf(2) = 1.0

        npr = 1
        pID(1) = 96
        pISO(1) = 2
        pf(1) = 1.0


        ref = 'N/A'

    end subroutine reaction0043_13c

    !**********************************************************************************************
    !(13C)H + H → H2 + (13C)
    subroutine reaction0044_13c(nh,p,t,dens,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf,ref)

        !==========================================================================================
        !       (13C)H + H → H2 + (13C)
        !==========================================================================================

        !Inputs
        integer, intent(in) :: nh                               
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Local
        integer :: rtype1,ns1,npr1,sID1(2),sISO1(2),pID1(2),pISO1(2)
        double precision :: sf1(2),pf1(2),rrates1(nh)
        character (len = 100) :: ref1

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
        character (len = 100) :: ref

        !Calling reaction r0043 to get reaction rate
        call reaction0044(nh,p,t,dens,rrates1,rtype1,ns1,sID1,sISO1,sf1,npr1,pID1,pISO1,pf1,ref1)

        rrates(:) = rrates1(:)*1.0

        rtype = 3

        ns = 3
        sID(1) = 93
        sISO(1) = 3
        sf(1) = 1.0
        sID(2) = 48
        sISO(2) = 0
        sf(2) = 1.0

        npr = 2
        pID(1) = 39
        pISO(1) = 0
        pf(1) = 1.0
        pID(2) = 46
        pISO(2) = 2
        pf(2) = 1.0


        ref = 'N/A'

    end subroutine reaction0044_13c

    !**********************************************************************************************



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
        double precision, intent(in) :: p(nh),t(nh),dens(nh)    

        !Local
        integer :: rtype1,ns1,npr1,sID1(2),sISO1(2),pID1(2),pISO1(2)
        double precision :: sf1(2),pf1(2),rrates1(nh)
        character (len = 100) :: ref1

        !Outputs
        double precision, intent(out) :: rrates(nh)             
        integer, intent(out) :: rtype
        integer, intent(out) :: ns,sID(2),sISO(2)
        integer, intent(out) :: npr,pID(2),pISO(2)
        double precision, intent(out) :: sf(2),pf(2)
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