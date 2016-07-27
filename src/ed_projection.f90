module ed_projection

    use mpi
    use dmft_grid, only: omega, nwloc
    use dmft_params, only: norb, mu, nw, nspin, nbath, nsite
    use dmft_green, only: cluster_hybridization_ftn

    public :: fit_h_imp_params

    double precision, allocatable :: &
        em_old(:,:)

    double complex, allocatable :: &
        D_cl_old(:,:,:), &
        G_cl_old(:,:,:), &
        G_loc_old(:,:,:)

    integer :: nx

    private
contains

    subroutine fit_h_imp_params( em, ek, vmk, D_cl, G_cl, G_loc )
        double precision, intent(inout) :: &
            em(norb,2), &
            ek(nbath,2), &
            vmk(norb,nbath,2)
        double complex, intent(in) :: &
            D_cl(nwloc,norb,nspin), &
            G_cl(nwloc,norb,nspin), &
            G_loc(nwloc,norb,nspin)

        double precision, allocatable :: x(:), df1(:), df2(:)
        integer :: ispin, iorb, ibath, ix, itr
        double precision :: tol, diff

        tol = 1.d-8

        nx = nspin*(norb + nbath + norb*nbath)

        allocate(D_cl_old(nwloc,norb,nspin), G_cl_old(nwloc,norb,nspin))
        allocate(G_loc_old(nwloc,norb,nspin), em_old(norb,nspin))
        allocate(x(nx),df1(nx),df2(nx))

        D_cl_old = D_cl
        G_cl_old = G_cl
        G_loc_old = G_loc
        em_old = em

        call emk_to_x( em, ek, vmk, x )

        ! ! test analytic derivative against numerical one
        ! print *, "em", em
        ! print *, "ek", ek
        ! print *, "vmk", vmk
        ! print *, "x", x
        ! print *, "func(x)", func(x)
        ! call dfunc(x,df1)
        ! call ndfunc(x,df2)
        ! print *, "derivative diff = ", sum(df1-df2)
        ! print *, "df1          df2"
        ! do ix=1,nx
        !     print *, df1(ix), df2(ix)
        ! enddo
        if (master) then
            write(*,"(1x,a,I6,a,E12.4)") &
                "Finding new em, ek, vmk ... "
        endif

        call FRPRMN( x, nx, tol, itr, diff )
        
        call x_to_emk( x, em, ek, vmk )

        if (master) then
            write(*,"(1x,a,I6,a,E12.4)") &
                "fitting procedure converged : itr=",itr,", diff=",diff
        endif

        ! for paramagnetic runs, set ek and vk spin independent
        if (nspin.eq.1) then
            ek(:,2) = ek(:,1)
            em(:,2) = em(:,1)
            vmk(:,:,2) = vmk(:,:,1)
        endif


        if (master) then
            write(*,*) 
            write(*,*) "New impurity Hamiltonian parameters"
            do ispin=1,nspin
                write(*,*) 
                write(*,*) "Spin ",ispin
                write(*,*) "Impurity/Bath Levels"
                do i=1,norb
                    write(*,"(1x,A,I2,F12.6)") "orb  ",i,em(i,ispin)
                enddo
                do i=1,nbath
                    write(*,"(1x,A,I2,F12.6)") "bath ",i,ek(i,ispin)
                enddo
                write(*,*)
                write(*,*) "Impurity/Bath Hybridization"
                write(*,"(7x)", advance="no")
                do i=1,norb
                    write(*,"(4x,A3,I1,4x)",advance="no") "orb",i
                enddo
                write(*,*)
                do i=1,nbath
                    write(*,"(1x,a4,I2)",advance="no") "bath",i
                    do j=1,norb
                        write(*,"(F12.6)",advance="no") vmk(j,i,ispin)
                    enddo
                    write(*,*)
                enddo
            enddo
            write(*,*)
        endif

        deallocate(D_cl_old, G_cl_old, G_loc_old, em_old)
        deallocate(x,df1,df2)
    end subroutine fit_h_imp_params

    subroutine emk_to_x( em, ek, vmk, x )
        double precision, intent(in) :: em(norb,2), ek(nbath,2), vmk(norb,nbath,2)
        double precision, intent(out) :: x(nx)
        integer :: iorb,ibath,ispin,ix
        ix=0
        do ispin=1,nspin
            do iorb=1,norb
                ix = ix+1
                x(ix) = em(iorb,ispin)
            enddo
            do ibath=1,nbath
                ix = ix+1
                x(ix) = ek(ibath,ispin)
            enddo
            do ibath=1,nbath
                do iorb=1,norb
                    ix = ix+1
                    x(ix) = vmk(iorb,ibath,ispin) 
                enddo
            enddo
        enddo
    end subroutine emk_to_x

    subroutine x_to_emk( x, em, ek, vmk )
        double precision, intent(in) :: x(nx)
        double precision, intent(out) :: em(norb,2), ek(nbath,2), vmk(norb,nbath,2)
        integer :: iorb,ibath,ispin,ix
        ix=0
        do ispin=1,nspin
            do iorb=1,norb
                ix = ix+1
                em(iorb,ispin) = x(ix)  
            enddo
            do ibath=1,nbath
                ix = ix+1
                ek(ibath,ispin) = x(ix) 
            enddo
            do ibath=1,nbath
                do iorb=1,norb
                    ix = ix+1
                    vmk(iorb,ibath,ispin) = x(ix)
                enddo
            enddo
        enddo
    end subroutine x_to_emk

    double precision function func(x)
        implicit none
        double precision :: x(nx)

        double precision :: func_loc, fr, fi
        integer :: iw, iorb, ibath, ispin

        double complex :: diff, D_cl_new(nwloc,norb,nspin)
        double precision :: &
            em(norb,2), &
            ek(nbath,2), &
            vmk(norb,nbath,2), &
            absdiff

        call x_to_emk( x, em, ek, vmk )
        call cluster_hybridization_ftn( ek, vmk, D_cl_new )

        func_loc = 0.0d0
        do ispin=1,nspin
            do iorb=1,norb
                do iw=1,nwloc
                    diff = em(iorb,ispin)-em_old(iorb,ispin) &
                     + D_cl_new(iw,iorb,ispin) - D_cl_old(iw,iorb,ispin) &
                     + 1.d0/G_loc_old(iw,iorb,ispin) - 1.d0/G_cl_old(iw,iorb,ispin)

                    absdiff = abs(diff)
                        
                    func_loc = func_loc + absdiff*absdiff*weight(iw)
                enddo
            enddo
        enddo

        call mpi_allreduce(func_loc,func,1,mpi_double_precision,&
                           mpi_sum,comm,mpierr)
        func = func/(norb*nspin)
    end function func

    double precision function weight(iw)
        integer :: iw
        weight = 1.d0/omega(iw)
    end function weight

    subroutine dfunc(x,df)
        implicit none
        double precision, intent(in) :: x(nx)
        double precision, intent(out) :: df(nx)
        
        double precision :: df_loc(nx)
        integer :: iw, iorb, ibath, ispin, ix

        double complex :: D_cl_new(nwloc,norb,nspin), &
                          diff(nwloc,norb,nspin), &
                          fac

        double precision :: &
            em(norb,2), &
            ek(nbath,2), &
            vmk(norb,nbath,2)

        call x_to_emk( x, em, ek, vmk )
        call cluster_hybridization_ftn( ek, vmk, D_cl_new )

        do ispin=1,nspin
            do iorb=1,norb
                do iw=1,nwloc
                    diff(iw,iorb,ispin) = em(iorb,ispin)-em_old(iorb,ispin) &
                     + D_cl_new(iw,iorb,ispin) - D_cl_old(iw,iorb,ispin) &
                     + 1.d0/G_loc_old(iw,iorb,ispin) - 1.d0/G_cl_old(iw,iorb,ispin)
                enddo
            enddo
        enddo


        do ispin=1,nspin
            ! df/dx w.r.t. em(iorb,ispin)
            do iorb=1,norb
                ix = (ispin-1)*(norb+nbath+norb*nbath) + iorb
                df_loc(ix) = 0.d0
                do iw=1,nwloc
                    df_loc(ix) = df_loc(ix) &
                                 + weight(iw)*2*real(diff(iw,iorb,ispin))
                enddo
            enddo
            ! df/dx w.r.t. ek(ibath,ispin)
            do ibath=1,nbath
                ix = (ispin-1)*(norb+nbath+norb*nbath) + norb + ibath
                
                df_loc(ix) = 0.d0

                do iorb=1,norb
                    do iw=1,nwloc
                        fac = vmk(iorb,ibath,ispin)**2
                        fac = fac/(cmplx(0.d0,-omega(iw))-ek(ibath,ispin))**2

                        df_loc(ix) = df_loc(ix) &
                            + weight(iw)*2*real(diff(iw,iorb,ispin)*fac)
                    enddo
                enddo
            enddo

            ! df/dx w.r.t. vmk(iorb,ibath,ispin)
            do ibath=1,nbath
                do iorb=1,norb
                    ix = (ispin-1)*(norb+nbath+norb*nbath) + norb + nbath &
                        + (ibath-1)*norb + iorb

                    df_loc(ix) = 0.d0

                    do iw=1,nwloc
                        fac = 2*vmk(iorb,ibath,ispin)
                        fac = fac/(cmplx(0.d0,-omega(iw))-ek(ibath,ispin))

                        df_loc(ix) = df_loc(ix) &
                            + weight(iw)*2*real(diff(iw,iorb,ispin)*fac)
                    enddo
                enddo
            enddo
        enddo

        call mpi_allreduce(df_loc,df,nx,mpi_double_precision,mpi_sum,comm,mpierr)
         df = df/(norb*nspin)
    end subroutine dfunc

    subroutine ndfunc(x,df)
        implicit none
        double precision, intent(in) :: x(nx)
        double precision, intent(out) :: df(nx)

        double precision, parameter :: h = 1.0d-6
        double precision :: f1, f2, xp(nx)
        integer :: i,j
        df = 0.0d0
        do i=1,nx
            xp = x
            xp(i) = xp(i)+h
            f1 = func(xp)
            xp(i) = xp(i)-2*h
            f2 = func(xp)

            df(i) = (f1-f2)/(2*h)
        enddo
    end subroutine ndfunc

    ! USES dfunc, func, linmin
    ! Given a starting point p that is a vector of length n, Fletcher-Reeves-Polak-Ribiere minimization
    ! is performed on a function func, using its gradient as calculated by a routine dfunc.
    ! The convergence tolerance on the function value is input as ftol. Returned quantities are p 
    ! (the location of the minimum), iter ( the # of iterations that were performed), and fret (the minimum value of function)
    ! The routine linmin is called to perform line minimizations
    ! parameter: nmax is the maximum anticipated value of n; itmax is the maximum allowd number of iterations;
    ! eps is a small number of rectify special case of converging to exactly zero function value.
    subroutine frprmn(p,n,ftol,iter,fret)
        implicit none

        integer:: iter, n, nmax, itmax,i
        double precision:: fret, ftol, p(n), eps
        parameter(nmax=60, itmax=20000, eps=1.d-10)

        integer its, j 
        double precision dgg, fp, gam, gg, g(nmax), h(nmax), xi(nmax),zero

        zero=0.D0

        fp = func(p)
        call dfunc(p,xi)

        do j = 1, n
            g(j) = -xi(j)
            h(j) = g(j)
            xi(j) = h(j)
        enddo

        do its = 1, itmax
            iter = its
            call linmin(p,xi,n,fret)
            if(2.D0*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+eps)) then
                return
            endif
            fp=fret
            call dfunc(p,xi)
            gg = 0.D0
            dgg = 0.D0

            do j = 1, n 
                gg = gg + g(j)*g(j)
                dgg = dgg+(xi(j)+g(j))*xi(j)
            enddo

            if(gg.eq.zero) return

            gam = dgg/gg
            do j =1, n
                g(j) = -xi(j)
                h(j) = g(j)+gam*h(j)
                xi(j)=h(j)
            enddo
        enddo

        stop "frprmn maximum iterations exceeded!"
        return
    end subroutine frprmn

    subroutine linmin(p,xi,n,fret)
        implicit none
        integer n,nmax  
        double precision fret,p(n),xi(n),tol
        parameter (nmax=60,tol=1.d-8)  
        !u    uses brent,f1dim,mnbrak  
        integer j,ncom  
        double precision ax,bx,fa,fb,fx,xmin,xx,pcom(nmax),xicom(nmax)
        do j=1,n  
            pcom(j)=p(j)  
            xicom(j)=xi(j)  
        enddo
        ax=0.d0  
        xx=1.d0
        call mnbrak(ax,xx,bx,fa,fx,fb,n,pcom,xicom)
        fret=brent(ax,xx,bx,tol,xmin,n,pcom,xicom)
        do j=1,n  
            xi(j)=xmin*xi(j)  
            p(j)=p(j)+xi(j)  
        enddo
        return  
    end subroutine linmin

    double precision function brent(ax,bx,cx,tol,xmin,n,pcom,xicom)
        implicit none
        integer itmax,n
        double precision ax,bx,cx,tol,xmin,cgold,zeps  
        external f  
        parameter (itmax=10000,cgold=.3819660,zeps=1.0d-10)  
        integer iter  
        double precision a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm  
        double precision :: pcom(n),xicom(n)

        a=min(ax,cx)  
        b=max(ax,cx)  
        v=bx  
        w=v  
        x=v  
        e=0.d0  
        fx=f1dim(x,n,pcom,xicom)
        fv=fx  
        fw=fx  
        do 11 iter=1,ITMAX  
            xm=0.5d0*(a+b)  
            tol1=tol*abs(x)+ZEPS  
            tol2=2.d0*tol1  
            if(abs(x-xm).le.(tol2-.5d0*(b-a))) goto 3  
            if(abs(e).gt.tol1) then  
                r=(x-w)*(fx-fv)  
                q=(x-v)*(fx-fw)  
                p=(x-v)*q-(x-w)*r  
                q=2.d0*(q-r)  
                if(q.gt.0.d0) p=-p  
                q=abs(q)  
                etemp=e  
                e=d  
                if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x))  goto 1    
                d=p/q  
                u=x+d  
                if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)  
                goto 2  
            endif  
1           if(x.ge.xm) then  
                e=a-x  
            else  
                e=b-x  
            endif  
            d=CGOLD*e  
2           if(abs(d).ge.tol1) then  
                u=x+d  
            else  
                u=x+sign(tol1,d)  
            endif  
            fu=f1dim(u,n,pcom,xicom)
            if(fu.le.fx) then  
                if(u.ge.x) then  
                    a=x  
                else  
                    b=x  
                endif  
                v=w  
                fv=fw  
                w=x  
                fw=fx  
                x=u  
                fx=fu  
            else  
                if(u.lt.x) then  
                    a=u  
                else  
                    b=u  
                endif  
                if(fu.le.fw .or. w.eq.x) then  
                    v=w  
                    fv=fw  
                    w=u  
                    fw=fu  
                else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then  
                    v=u  
                    fv=fu  
                endif  
            endif  
11      continue  
        stop 'brent exceed maximum iterations'  

3       xmin=x  
        brent=fx  
        return  
    end function brent

    double precision function f1dim(x,n,pcom,xicom)

        implicit none
        integer n
        double precision x  
        integer j
        double precision pcom(n),xicom(n),xt(n)

        do j=1,n
            xt(j)=pcom(j)+x*xicom(j)  
        enddo

        f1dim=func(xt)
        return  
    end function f1dim

    subroutine mnbrak(ax,bx,cx,fa,fb,fc,n,pcom,xicom)
        implicit none

        integer :: n
        double precision ax,bx,cx,fa,fb,fc,GOLD,GLIMIT,TINY  
        PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.D-20)  
        double precision dum,fu,q,r,u,ulim,pcom(n),xicom(n)

        fa=f1dim(ax,n,pcom,xicom)
        fb=f1dim(bx,n,pcom,xicom)  

        if(fb.gt.fa)then  
            dum=ax  
            ax=bx  
            bx=dum  
            dum=fb  
            fb=fa  
            fa=dum  
        endif  
        cx=bx+GOLD*(bx-ax)  
        fc=f1dim(cx,n,pcom,xicom)
1       if(fb.ge.fc)then  
            r=(bx-ax)*(fb-fc)  
            q=(bx-cx)*(fb-fa)  
            u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))  
            ulim=bx+GLIMIT*(cx-bx)  
            if((bx-u)*(u-cx).gt.0.d0)then  
                fu=f1dim(u,n,pcom,xicom)
                if(fu.lt.fc)then  
                    ax=bx  
                    fa=fb  
                    bx=u  
                    fb=fu  
                    return  
                else if(fu.gt.fb)then  
                    cx=u  
                    fc=fu  
                    return  
                endif  
                u=cx+GOLD*(cx-bx)  
                fu=f1dim(u,n,pcom,xicom)
            else if((cx-u)*(u-ulim).gt.0.)then  
                fu=f1dim(u,n,pcom,xicom)
                if(fu.lt.fc)then  
                    bx=cx  
                    cx=u  
                    u=cx+GOLD*(cx-bx)  
                    fb=fc  
                    fc=fu  
                    fu=f1dim(u,n,pcom,xicom)
                endif  
            else if((u-ulim)*(ulim-cx).ge.0.)then  
                u=ulim  
                fu=f1dim(u,n,pcom,xicom)
            else  
                u=cx+GOLD*(cx-bx)  
                fu=f1dim(u,n,pcom,xicom)
            endif  
            ax=bx  
            bx=cx  
            cx=u  
            fa=fb  
            fb=fc  
            fc=fu  
            goto 1  
        endif  
        return  
    end subroutine
end module ed_projection
