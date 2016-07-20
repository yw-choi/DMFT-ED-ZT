module ed_projection

    use dmft_grid, only: omega, nwloc
    use ed_params, only: nbath
    use dmft_params, only: norb, mu, nw, nspin
    use mpi

    public :: project_to_impurity_model

    integer :: nx, iorb, np

    double complex, allocatable :: G0(:)

    private
contains

    ! obtain ek,vk by minimizing the difference between the cluster Weiss field and G0
    subroutine project_to_impurity_model(G0in,ek,vk)
        implicit none
        double complex, intent(in) :: G0in(nwloc,norb,nspin)
        double precision, intent(inout) :: ek(norb+nbath,2), vk(norb,nbath,2)

        ! local variables
        integer :: itr, i,j, ispin
        double precision, parameter :: tol = 1.0D-7 ! minimization tolerance 
        double precision :: x(1+2*nbath/norb), diff

        double precision :: df1(1+2*nbath/norb), df2(1+2*nbath/norb)

        np = nbath/norb
        ! # of fitting parameters for each orbital
        nx = 1 + 2*np

        allocate(G0(nwloc))

        if (master) then
            write(*,*)
            write(*,*) "Projecting to the impurity model..."
        endif

        do ispin=1,nspin
            do i=1,norb
                ! these following two variables are module-wide variables
                ! used in other subroutines
                iorb = i
                G0 = G0in(:,iorb,ispin)

                ! collect fitting parameters
                x(1) = ek(iorb,ispin)
                do j=1,np
                    x(1+j) = ek(norb*j+iorb,ispin)
                    x(1+np+j) = vk(iorb,norb*(j-1)+iorb,ispin)
                enddo

                ! find x that minimizing the differences b/w G0 and G0cl
                call FRPRMN(x,nx,tol,itr,diff)

                ! ! test analytic derivative against numerical one
                ! call dfunc(x,df1)
                ! call ndfunc(x,df2)
                ! print *, "derivative diff = ", sum(df1-df2)

                ! recast x to ek,vk
                ek(iorb,ispin) = x(1)

                do j=1,np
                    ek(norb*j+iorb,ispin) = x(1+j)
                    vk(iorb,norb*(j-1)+iorb,ispin) = x(1+np+j)
                enddo

                if (master) then
                    write(*,"(1x,a,I1,a,I1,a,I6,a,E12.4)") &
                        "(iorb,ispin)=(",iorb,",",ispin,") converged : itr=",&
                        itr,", diff=",diff
                endif

            enddo
        enddo

        ! for paramagnetic runs, set ek and vk spin independent
        if (nspin.eq.1) then
            ek(:,2) = ek(:,1)
            vk(:,:,2) = vk(:,:,1)
        endif

        if (master) then
            write(*,*) 
            write(*,*) "Projected impurity/bath Levels"
            do ispin=1,nspin
                write(*,*) 
                write(*,*) "Spin ",ispin
                do i=1,norb
                    write(*,"(1x,A,I2,F12.6)") "orb ",i,ek(i,ispin)
                enddo
                do i=norb+1,norb+nbath
                    write(*,"(1x,A,I2,F12.6)") "bath ",(i-norb),ek(i,ispin)
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
                        write(*,"(F12.6)",advance="no") vk(j,i,ispin)
                    enddo
                    write(*,*)
                enddo
            enddo
            write(*,*)
        endif

        deallocate(G0)
    end subroutine project_to_impurity_model

    ! cluster hybridization function evaluation 
    double complex function delta_cl(iw,x)
        integer :: iw
        double precision :: x(nx)

        integer :: i, xidx

        delta_cl = cmplx(0.0d0,0.0d0)

        do i=1,np
            delta_cl = delta_cl + x(1+np+i)*x(1+np+i)/(cmplx(0.0d0,omega(iw))-x(1+i))
        enddo

    end function delta_cl

    double precision function func(x)
        implicit none
        double precision :: x(nx)
        double precision :: func_loc, diff
        double complex :: gf

        integer :: iw, ibath, xidx

        func_loc = 0.0d0
        do iw=1,nwloc
            gf = g0cl(iw,x)
            diff = abs(gf - G0(iw))

            func_loc = func_loc + diff*diff*weight(iw)
        enddo

        call mpi_allreduce(func_loc,func,1,mpi_double_precision,mpi_sum,comm,mpierr)
    end function func

    double precision function weight(iw)
        integer :: iw
        weight = 1 !/omega(iw)
    end function

    double complex function g0cl(iw,x)
        integer :: iw
        double precision :: x(nx)

        g0cl = cmplx(0.d0,omega(iw)) + mu - x(1) - delta_cl(iw,x)

        g0cl = 1.0d0/g0cl
    end function g0cl

    subroutine ndfunc(x,df)
        implicit none
        double precision, intent(in) :: x(nx)
        double precision, intent(out) :: df(nx)

        double precision, parameter :: h = 1.0d-8
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

    ! numerical derivative
    ! subroutine dfunc(x,df)
    !     implicit none
    !     double precision, intent(in) :: x(nx)
    !     double precision, intent(out) :: df(nx)

    !     double precision, parameter :: h = 1.0d-6
    !     double precision :: f1, f2, xp(nx)
    !     integer :: i,j
    !     df = 0.0d0
    !     do i=1,nx
    !         xp = x
    !         xp(i) = xp(i)+h
    !         f1 = func(xp)
    !         xp(i) = xp(i)-2*h
    !         f2 = func(xp)

    !         df(i) = (f1-f2)/(2*h)
    !     enddo
    ! end subroutine dfunc

    subroutine dfunc(x,df)
        implicit none
        double precision, intent(in) :: x(nx)
        double precision, intent(out) :: df(nx)
        
        double precision :: df_loc(nx), fac
        double complex :: gf, diffconjg, fac2, fac3
        integer :: iw, i

        df_loc = 0.0d0

        do iw=1,nwloc
            gf = g0cl(iw,x)
            diffconjg = conjg(gf-g0(iw))

            fac = 2*weight(iw)

            fac2 = diffconjg*gf*gf

            df_loc(1) = df_loc(1)+fac*real(fac2)

            do i=1,np
                fac3 = x(1+np+i)*x(1+np+i)/(cmplx(0.d0,omega(iw))-x(1+i))**2
                df_loc(1+i) = df_loc(1+i)+fac*real(fac2*fac3)

                fac3 = 2*x(1+np+i)/(cmplx(0.d0,omega(iw))-x(1+i))
                df_loc(1+np+i) = df_loc(1+np+i)+fac*real(fac2*fac3)
            enddo
        enddo

        call mpi_allreduce(df_loc,df,nx,mpi_double_precision,mpi_sum,comm,mpierr)
    end subroutine dfunc

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

        ! if (master) then
        !     write(*,*) "initial value and differences"
        !     write(*,*) "x                       df/dx"
        !     do j=1,nx
        !        write(*,*) p(j), xi(j) 
        !     enddo
        ! endif

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
