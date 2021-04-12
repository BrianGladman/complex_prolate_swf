program cprodrv
    use param
    use complex_prolate_swf
    
    implicit none
    integer      mmin, minc, mnum, lnum, ioprad, iopang, iopnorm, minacc, &
                 ngau, ioparg, narg, maxd, maxdr, maxint, maxj, maxlp, maxm, & 
                 maxmp, maxn, maxp, maxpdr, maxq, maxt, neta, jnenmax, ndec, &
                 nex, kindd, kindq, im, i, j, l, m, maxc, lnump

    real(knd)    arg1, c, darg, x1
    complex(knd) cc

    real (knd), dimension(:), allocatable ::        eta
    complex(knd), dimension(:), allocatable ::      r1c, r1dc, r2c, r2dc
    integer, dimension(:), allocatable ::           ir1e, ir1de, ir2e, ir2de, naccr
    complex(knd), dimension (:,:), allocatable ::   s1c, s1dc
    integer, dimension(:,:), allocatable ::         is1e, is1de, naccs, naccds

    kindd =  8
    kindq = 16

    open(1, file='cprofcn.dat')
    open(20, file='fort.20')
    open(30, file='fort.30')
    open(40, file='fort.20')
    open(50, file='fort.30')
    open(60, file='fort.20')

    read(1,*) mmin, minc, mnum, lnum
    read(1,*) ioprad, iopang, iopnorm
    read(1,*) cc, x1
    if(iopang /= 0) read(1,*) ioparg, arg1, darg, narg
        
    allocate (eta(narg), r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate (ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    allocate (s1c(lnum, narg), s1dc(lnum, narg))
    allocate (is1e(lnum, narg), is1de(lnum, narg), naccs(lnum, narg), naccds(lnum, narg))

    kindd = 8
    kindq = 16
    if(knd == kindd) minacc = 8
    if(knd == kindq .and. aimag(cc) <= 20.0e0_knd) minacc = 15
    if(knd == kindq .and. aimag(cc) > 20.0e0_knd) minacc = 8
    c = abs(cc)
    maxc = 2 * int((x1 + 1.0e0_knd) * (abs(real(cc)) + abs(aimag(cc)))) + 25
    maxm = mmin + minc * (mnum - 1)
    maxint = lnum + 3 * ndec + int(c) + 5
    maxj = maxint + maxm
    maxp = maxint
    maxn = maxp + maxm
    maxn = max(maxn, maxc)
    maxpdr = 4 * ndec + 5
    neta = 993
    ngau = 200
    if(ioprad == 2) then
        lnump = max(lnum + maxm, 1000)
        if(x1 >= 0.00065e0_knd) maxn = 2 * (lnump * (-18.5e0_knd - 20.0e0_knd * log10(x1)) & 
                                        + 5 * ndec + 4 * maxm + c + 5000) + maxm + 5
        if(x1 > 0.08e0_knd) maxn = 2 * (lnump * (0.5e0_knd - 3.0e0_knd * log10(x1)) &
                                        + 5 * ndec + 4 * maxm + c + 1000) + maxm + 5
        if(x1 > 1.0e0_knd) maxn = 2 * (lnump * 0.5e0_knd + 5 * ndec + 4 * maxm + c + 500) + maxm + 5
        if(x1 <= 0.5e0_knd) maxpdr = maxpdr + int(2.e0_knd * c + 100.0e0_knd * x1) + 400
        maxn = max(maxn, maxc)
        maxp = max(maxn, maxp)
        if(x1 < 1.0e-3_knd) ngau = 200 - 50 * int(log10(x1) - 1.0e-30_knd)
        if(x1 < 1.0e-10_knd) ngau = 250 - 50 * int(log10(x1) - 1.0e-30_knd)
        if(x1 < 1.0e-11_knd) ngau = 1200
        if(x1 < 1.0e-12_knd) ngau = 2400
    end if
    maxq = maxint + maxm + maxm
    maxd = maxp / 2 + 1
    maxdr = maxpdr / 2 + 1
    maxp = max(maxp, maxpdr) + 5
    maxlp = lnum + maxm + 5
    maxmp = maxm + 5
    maxt = 1
    jnenmax = 10
    if(iopang /= 0) maxt = narg
    
    if (iopang /= 0) then
        do j = 1, narg  
            eta(j) = arg1 + (j - 1) * darg
        end do
    end if

    do im = 1, mnum
        m = mmin + (im - 1) * minc

        call cprofcn(cc, m, lnum, ioprad, x1, iopang, iopnorm, narg, eta, &
                       r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                       s1c, is1e, s1dc, is1de, naccs, naccds)

        if (ioprad /= 0) then

            if(knd == kindd .and. ioprad /= 0) write(20, 120) x1 + 1.0e0_knd, cc, m
120         format(1x,'x = 'e23.14,'; c = ',e23.14, e23.14,'; m = ',i5)
            if(knd == kindq .and. ioprad /= 0) write(20, 130) x1 + 1.0e0_knd, cc, m
130         format(1x,'x = 'e39.30,'; c = ',e39.30, e39.30,'; m = ',i5)
            
            do i = 1, lnum
                l = m + i - 1
                
                if(ioprad == 1) write(20, 186) l, r1c(i), ir1e(i), r1dc(i), ir1de(i), naccr(i)
186             format(1x, i5, 2x, 2(f17.14, 1x, f17.14, i5, 2x), i4)
                
                if(ioprad == 2) write(20, 690) l, r1c(i), ir1e(i), r1dc(i), ir1de(i), r2c(i), ir2e(i), r2dc(i), ir2de(i), naccr(i)
690             format(1x, i5, 2x, 2(f17.14, 1x, f17.14, i5, 2x),/,8x, 2(f17.14, 1x, f17.14, i5, 2x), i2, '  ')
                
            end do
        end if
            
        if (iopang /= 0) then
            
            if(knd == kindd .and. iopang /= 0) write(30, 65) cc, m
65          format(1x,'c = ',e23.14, e23.14,'; m = ',i5)
            if(knd == kindq .and. iopang /= 0) write(30, 70) cc, m
70          format(1x,'c = ',e39.30, e39.30,'; m = ',i5)
                
            do i = 1, lnum
                l = m + i - 1

                write(30, "(1x,i6)") l

                do j = 1, narg

                    if(iopang == 1) write(30, 750) eta(j), s1c(i,j), is1e(i,j), naccs(i,j)
                    if(iopang == 2) write(30, 760) eta(j), s1c(i,j), is1e(i,j), s1dc(i,j), is1de(i,j), naccs(i,j), naccds(i,j)
750                 format(1x, f19.14, 2x, f17.14, 1x, f17.14, 2x, i5, 2x, i2)
760                 format(1x, f19.14, 2x, f17.14, 1x, f17.14, 2x, i5, 2x, f17.14, 1x, f17.14, 2x, i5, 2x, i2, ', ', i2)
                
                end do
            end do
        end if
    end do
    
    deallocate (is1e, is1de, naccs, naccds)
    deallocate (s1c, s1dc)
    deallocate (ir1e, ir1de, ir2e, ir2de, naccr)
    deallocate (eta, r1c, r1dc, r2c, r2dc)
    close(60)
    close(50)
    close(40)
    close(30)
    close(20)
    close(1)

end program cprodrv
