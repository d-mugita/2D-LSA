program NELFA
    implicit none
!================ Definition
!---- physical constant
    double precision, parameter :: Pi = 3.141592653589793d0

!---- system variables
    double precision, parameter :: lx = 64.d0 ! box size
    double precision, parameter :: ly = lx
    integer, parameter          :: num_part = 4096 ! total particle number
    double precision, parameter :: pack_frac = 0.720d0 ! packing fraction
    double precision, parameter :: rad_ratio = 1.4d0 ! radius ratio
    double precision, parameter :: mol_ratio = 1.d0/3 ! molar ratio
    integer, parameter          :: num_part_L = int(num_part * mol_ratio) ! large particle number
    integer, parameter          :: num_part_S = num_part - num_part_L ! small particle number
    double precision, parameter :: rad_S = dsqrt(lx*ly*pack_frac/((dble(num_part_S)+num_part_L*rad_ratio**2)*Pi)) ! small particle radius
    double precision, parameter :: rad_L = rad_S * rad_ratio ! large particle radius
    double precision, dimension(2*num_part) :: pos ! particle position
    double precision, dimension(num_part)   :: rad_list ! radius list
    double precision, dimension(num_part)   :: mass_list ! mass list

!---- grid mapping
    integer, parameter                     :: gd_num_x = int(lx/(dsqrt(2.d0)*rad_S)+1)
    integer, parameter                     :: gd_num_y = int(ly/(dsqrt(2.d0)*rad_S)+1)
    integer, dimension(gd_num_x,gd_num_y)  :: gd_map
    integer, dimension(2,num_part)         :: gd_coord
    integer, parameter                     :: mask_num = 120 ! grid length 5: 120, 9: 292
    integer, dimension(2,mask_num)         :: mask
    integer, dimension(mask_num,num_part)  :: gd_neigh_list

!---- NELF-A
    integer, parameter :: cand_num_max = 40 ! mono: 20, bin: 40
    integer, parameter :: FP_num_max = cand_num_max ! FP: Free volume construction Particles
    double precision, dimension(num_part) :: FV, FS ! FV: Free volume, FS: Free Surface area

!---- other
    integer(8) :: i

!================ Initialization
    open(10, file='./bin_4096_720.dat')
    do i = 1, num_part
        read(10,*) pos(i), pos(i+num_part)
    end do
    close(10)  
    call make_rad_list(ly,num_part,num_part_L,rad_S,rad_L,rad_list)
    mass_list = 1.0D0
    call gridmask5(mask)

!================ Main
    call gridmapping(lx,ly,num_part,pos,gd_num_x,gd_num_y,gd_map,gd_coord)
    call finding_neighbors(num_part,gd_num_x,gd_num_y,gd_map,gd_coord,mask_num,mask,gd_neigh_list)

    call NELFA_cutoff(Pi,lx,ly,num_part,pos,rad_list,mask_num,gd_neigh_list,&
                      cand_num_max,FP_num_max,FV,FS)

    open(20, file='./FV.dat')
    open(21, file='./FS.dat')
    do i = 1, num_part
        write(20,*) FV(i) 
        write(21,*) FS(i)
    end do
    close(20)
    close(21)
end

!================ Subroutine
subroutine make_rad_list(ly,num_part,num_part_L,rad_S,rad_L,rad_list)
    implicit none
    double precision, intent(in) :: ly
    integer, intent(in) :: num_part
    integer, intent(in) :: num_part_L
    double precision, intent(in) :: rad_S
    double precision, intent(in) :: rad_L
    double precision, dimension(num_part), intent(out) :: rad_list
    integer :: num_L_max
    integer :: i, k

    num_L_max = 0
    do i = 1, num_part-1, 2
        k = i/int(ly)+1
        if ( mod(k,2) == 1 ) then
            rad_list(i) = rad_S
            if ( num_L_max <= num_part_L ) then
                rad_list(i+1) = rad_L
                num_L_max = num_L_max+1
            else
                rad_list(i+1) = rad_S
            end if
        else
            if ( num_L_max <= num_part_L ) then
                rad_list(i) = rad_L
                num_L_max = num_L_max+1
            else
                rad_list(i) = rad_S
            end if
            rad_list(i+1) = rad_S
        end if
    end do
    return
end subroutine

subroutine gridmask2(mask)
    implicit none
    integer, dimension(2,24), intent(out) :: mask
    integer, dimension(24) :: unit_mask_x, unit_mask_y
    integer :: i

    data unit_mask_x/ -2,-2,-2,-2,-2,&
                      -1,-1,-1,-1,-1,&
                       0, 0,    0, 0,&
                       1, 1, 1, 1, 1,&
                       2, 2, 2, 2, 2/
    data unit_mask_y/ -2,-1, 0, 1, 2,&
                      -2,-1, 0, 1, 2,&
                      -2,-1,    1, 2,&
                      -2,-1, 0, 1, 2,&
                      -2,-1, 0, 1, 2/
    do i = 1, 24
        mask(1,i) = unit_mask_x(i)
        mask(2,i) = unit_mask_y(i)
    end do
    return
end subroutine

subroutine gridmask3(mask)
    implicit none
    integer, dimension(2,48), intent(out) :: mask
    integer, dimension(48) :: unit_mask_x, unit_mask_y
    integer :: i

    data unit_mask_x/ -3,-3,-3,-3,-3,-3,-3,&
                      -2,-2,-2,-2,-2,-2,-2,&
                      -1,-1,-1,-1,-1,-1,-1,&
                       0, 0, 0,    0, 0, 0,&
                       1, 1, 1, 1, 1, 1, 1,&
                       2, 2, 2, 2, 2, 2, 2,&
                       3, 3, 3, 3, 3, 3, 3/
    data unit_mask_y/ -3,-2,-1, 0, 1, 2, 3,&
                      -3,-2,-1, 0, 1, 2, 3,&
                      -3,-2,-1, 0, 1, 2, 3,&
                      -3,-2,-1,    1, 2, 3,&
                      -3,-2,-1, 0, 1, 2, 3,&
                      -3,-2,-1, 0, 1, 2, 3,&
                      -3,-2,-1, 0, 1, 2, 3/
    do i = 1, 48
        mask(1,i) = unit_mask_x(i)
        mask(2,i) = unit_mask_y(i)
    end do
    return
end subroutine

subroutine gridmask5(mask)
    implicit none
    integer, dimension(2,120), intent(out) :: mask
    integer, dimension(120) :: unit_mask_x, unit_mask_y
    integer :: i

    data unit_mask_x/ -5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,&
                      -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,&
                      -3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,&
                      -2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,&
                      -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,&
                       0, 0, 0, 0, 0,    0, 0, 0, 0, 0,&
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,&
                       2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,&
                       3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,&
                       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,&
                       5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5/
    data unit_mask_y/ -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,&
                      -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,&
                      -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,&
                      -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,&
                      -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,&
                      -5,-4,-3,-2,-1,    1, 2, 3, 4, 5,&
                      -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,&
                      -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,&
                      -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,&
                      -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,&
                      -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5/
    do i = 1, 120
        mask(1,i) = unit_mask_x(i)
        mask(2,i) = unit_mask_y(i)
    end do
    return
end subroutine

subroutine gridmask9(mask)
    implicit none
    integer, dimension(2,292), intent(out) :: mask
    integer, dimension(292) :: unit_mask_x, unit_mask_y
    integer :: i

    data unit_mask_x/         -9,-9,-9,-9,-9,-9,-9,&
                        -8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,&
                     -7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,&
                  -6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,&
               -5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,&
               -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,&
            -3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,&
            -2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,&
            -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,&
             0, 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0, 0,&
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,&
             2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,&
             3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,&
                4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,&
                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,&
                   6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,&
                      7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,&
                         8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,&
                               9, 9, 9, 9, 9, 9, 9/
    data unit_mask_y/         -3,-2,-1, 0, 1, 2, 3,&
                        -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,&
                     -6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6,&
                  -7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7,&
               -8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,&
               -8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,&
            -9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,&
            -9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,&
            -9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,&
            -9,-8,-7,-6,-5,-4,-3,-2,-1,    1, 2, 3, 4, 5, 6, 7, 8, 9,&
            -9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,&
            -9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,&
            -9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,&
               -8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,&
               -8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,&
                  -7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7,&
                     -6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6,&
                        -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,&
                              -3,-2,-1, 0, 1, 2, 3/
    do i = 1, 292
        mask(1,i) = unit_mask_x(i)
        mask(2,i) = unit_mask_y(i)
    end do
    return
end subroutine

subroutine gridmapping(lx,ly,num_part,pos,gd_num_x,gd_num_y,gd_map,gd_coord)
    implicit none
    double precision, intent(in) :: lx, ly
    integer, intent(in) :: num_part
    double precision, dimension(2*num_part), intent(in) :: pos
    integer, intent(in) :: gd_num_x, gd_num_y
    integer, dimension(gd_num_x,gd_num_y), intent(out) :: gd_map
    integer, dimension(2,num_part), intent(out) :: gd_coord
    integer :: i
    integer :: gd_x, gd_y
    double precision :: ilx, ily 

    gd_map = 0
    ilx = 1.d0/lx * gd_num_x
    ily = 1.d0/ly * gd_num_y
    do i = 1, num_part
        gd_x = int(ilx * pos(i)) + 1
        gd_y = int(ily * pos(i+num_part)) + 1
        gd_map(gd_x,gd_y) = i
        gd_coord(1,i) = gd_x
        gd_coord(2,i) = gd_y
    end do
    return
end subroutine

subroutine finding_neighbors(num_part,gd_num_x,gd_num_y,gd_map,gd_coord,mask_num,mask,gd_neigh_list)
    implicit none
    integer, intent(in) :: num_part
    integer, intent(in) :: gd_num_x, gd_num_y
    integer, dimension(gd_num_x,gd_num_y), intent(in) :: gd_map
    integer, dimension(2,num_part), intent(in) :: gd_coord
    integer, intent(in) :: mask_num
    integer, dimension(2,mask_num), intent(in) :: mask
    integer, dimension(mask_num,num_part), intent(out) :: gd_neigh_list
    integer :: mask_num2
    integer :: gd_x, gd_y
    integer :: i, j, k, m

    mask_num2 = int(mask_num/2)
    do i = 1, num_part
        m = 1
        do j = 1, mask_num2
            gd_x = gd_coord(1,i) + mask(1,j)
            gd_y = gd_coord(2,i) + mask(2,j)
            if ( gd_x > gd_num_x ) then
                gd_x = gd_x - gd_num_x
            else if ( gd_x <= 0 ) then
                gd_x = gd_x + gd_num_x
            end if
            if ( gd_y > gd_num_y ) then
                gd_y = gd_y - gd_num_y
            else if ( gd_y <= 0 ) then
                gd_y = gd_y + gd_num_y
            end if
            k = gd_map(gd_x,gd_y)
            if ( k /= 0 ) then
                gd_neigh_list(m,i) = k
                m = m + 1
            end if
        end do
        gd_neigh_list(mask_num-1,i) = m - 1
        do j = mask_num2+1, mask_num
            gd_x = gd_coord(1,i) + mask(1,j)
            gd_y = gd_coord(2,i) + mask(2,j)
            if ( gd_x > gd_num_x ) then
                gd_x = gd_x - gd_num_x
            else if ( gd_x <= 0 ) then
                gd_x = gd_x + gd_num_x
            end if
            if( gd_y > gd_num_y ) then
                gd_y = gd_y - gd_num_y
            else if ( gd_y <= 0 ) then
                gd_y = gd_y + gd_num_y
            end if
            k = gd_map(gd_x,gd_y)
            if ( k /= 0 ) then
                gd_neigh_list(m,i) = k
                m = m + 1
            end if
        end do
        gd_neigh_list(mask_num,i) = m - 1
    end do
    return
end subroutine

subroutine sort_integer_nochange(N,M,A,B)
    implicit none
    integer, intent(in) :: N ! array size
    integer, intent(in) :: M ! sorting part
    double precision, dimension(N), intent(in) :: A ! standard array (no change)
    integer, dimension(N), intent(inout) :: B ! array
    double precision, dimension(N) :: dummy
    double precision :: tmp1
    integer :: tmp2
    integer :: i, j, k

    dummy = A
    do i = 2, M
        tmp1 = dummy(i)
        tmp2 = B(i)
        do j = i-1, 1, -1
            if ( dummy(j) > tmp1 ) then
                dummy(j+1) = dummy(j)
                B(j+1) = B(j)
                k = j
            else
                k = j + 1
                exit
            end if
        end do
        dummy(k) = tmp1
        B(k) = tmp2
    end do
    return
end subroutine

subroutine sort_double_nochange(N,M,A,B)
    implicit none
    integer, intent(in) :: N ! array size
    integer, intent(in) :: M ! sorting part
    double precision, dimension(N), intent(in) :: A ! standard array (no change)
    double precision, dimension(N), intent(inout) :: B ! array
    double precision, dimension(N) :: dummy
    double precision :: tmp1
    double precision :: tmp2
    integer :: i, j, k

    dummy = A
    do i = 2, M
        tmp1 = dummy(i)
        tmp2 = B(i)
        do j = i-1, 1, -1
            if ( dummy(j) > tmp1 ) then
                dummy(j+1) = dummy(j)
                B(j+1) = B(j)
                k = j
            else
                k = j + 1
                exit
            end if
        end do
        dummy(k) = tmp1
        B(k) = tmp2
    end do
    return
end subroutine
  
subroutine NELFA_cutoff(Pi,lx,ly,num_part,pos,rad_list,mask_num,gd_neigh_list,&
                        cand_num_max,FP_num_max,FV,FS) 
    implicit none
    double precision, intent(in) :: Pi
    double precision, intent(in) :: lx, ly
    integer, intent(in) :: num_part
    double precision, dimension(2*num_part), intent(in) :: pos
    double precision, dimension(num_part), intent(in) :: rad_list
    integer, intent(in) :: mask_num
    integer, dimension(mask_num,num_part), intent(in) :: gd_neigh_list
    integer, intent(in) :: cand_num_max
    integer, intent(in) :: FP_num_max
    double precision, dimension(num_part), intent(out) :: FV, FS

    ! IS: InterSection
    ! FP, FI: Free volume construction Particles and Intersections
    double precision :: xi, yi
    double precision :: radi
    double precision :: radi2
    double precision :: min_dis
    integer :: FP1
    double precision :: dis_x, dis_y
    double precision :: dis_j
    double precision :: x1_reli, y1_reli
    integer, dimension(2*cand_num_max) :: sur_part
    double precision, dimension(4*cand_num_max) :: relpos_sur_IS
    integer :: len_sur_part
    double precision, dimension(2*cand_num_max) :: dis2_IS
    integer, dimension(1) :: ind_FP2
    integer :: FP2
    double precision :: x2_reli, y2_reli
    double precision :: x_IS, y_IS
    double precision :: vec1x, vec1y, vec2x, vec2y
    integer, dimension(FP_num_max) :: FP
    double precision, dimension(2*FP_num_max) :: pos_FP
    double precision, dimension(2*FP_num_max) :: pos_FI
    integer :: start
    integer :: con1, con2
    double precision :: x_con2, y_con2
    double precision :: dis12
    integer :: ind12
    integer :: ind_next
    integer :: next_part
    integer :: num_FP
    double precision :: poly_area, seg_area, seg_len
    integer :: j_next
    double precision :: angle
    integer i, j, k, m

    do i = 1, num_part
        xi = pos(i)
        yi = pos(i+num_part)
        radi = rad_list(i)
        radi2 = 2.d0*radi
        ! calc. FP1: 1st FP (= The nearest particle)
        min_dis = 100.d0*radi
        FP1 = 0
        do m = 1, gd_neigh_list(mask_num,i)
            j = gd_neigh_list(m,i)
            if ( j /= 0 ) then
                dis_x = pos(j) - pos(i)
                dis_y = pos(j+num_part) - pos(i+num_part)
                if ( dis_x >= lx/2 ) then
                    dis_x = dis_x - lx
                else if ( dis_x < -lx/2 ) then
                    dis_x = dis_x + lx
                end if
                if ( dis_y >= ly/2 ) then
                    dis_y = dis_y - ly
                else if ( dis_y < -ly/2 ) then
                    dis_y = dis_y + ly
                end if
                dis_j = dsqrt(dis_x**2 + dis_y**2) - rad_list(j)
                if ( dis_j < min_dis ) then
                    min_dis = dis_j
                    FP1 = j
                    x1_reli = dis_x
                    y1_reli = dis_y
                end if
            end if
        end do
        ! calc. the nearest IS and FP2
        call calc_sur_IS(FP1,i,lx,ly,num_part,pos,rad_list,mask_num,gd_neigh_list,&
                         cand_num_max,radi,radi2,sur_part,relpos_sur_IS,len_sur_part)
        do j = 1, len_sur_part
            dis2_IS(j) = (relpos_sur_IS(j) + x1_reli)**2 + &
                         (relpos_sur_IS(j+2*cand_num_max) + y1_reli)**2
        end do
        ind_FP2 = minloc(dis2_IS(:len_sur_part))
        FP2 = sur_part(ind_FP2(1))
        x2_reli = pos(FP2) - xi
        y2_reli = pos(FP2+num_part) - yi
        if ( x2_reli >= lx/2 ) then
            x2_reli = x2_reli - lx
        else if ( x2_reli < -lx/2 ) then
            x2_reli = x2_reli + lx
        end if
        if ( y2_reli >= ly/2 ) then
            y2_reli = y2_reli - ly
        else if ( y2_reli < -ly/2 ) then
            y2_reli = y2_reli + ly
        end if
        x_IS = relpos_sur_IS(ind_FP2(1)) + x1_reli
        y_IS = relpos_sur_IS(ind_FP2(1)+2*cand_num_max) + y1_reli
        ! decide con1, con2 (con: confirmed particle)
        vec1x = x1_reli - x_IS ! vec1: vector from IS to FP1
        vec1y = y1_reli - y_IS
        vec2x = x2_reli - x_IS ! vec1: vector from IS to FP2
        vec2y = y2_reli - y_IS
        if ( vec1x * vec2y - vec1y * vec2x > 0 ) then
            FP(1) = FP1
            pos_FP(1) = x1_reli
            pos_FP(1+FP_num_max) = y1_reli
            FP(2) = FP2
            pos_FP(2) = x2_reli
            pos_FP(2+FP_num_max) = y2_reli
            pos_FI(1) = x_IS
            pos_FI(1+FP_num_max) = y_IS
            start = FP1
            con1 = FP1
            con2 = FP2
            x_con2 = x2_reli
            y_con2 = y2_reli
            ! calc. surrounding IS and the index of intersection
            call calc_sur_IS(con2,i,lx,ly,num_part,pos,rad_list,mask_num,gd_neigh_list,&
                             cand_num_max,radi,radi2,sur_part,relpos_sur_IS,len_sur_part)
            dis12 = -1.d0
            do k = 1, len_sur_part
                if ( sur_part(k) == con1 ) then
                    if ( dis12 == -1.d0 ) then
                        dis12 = (x_IS - relpos_sur_IS(k) - x_con2)**2 + &
                                (y_IS - relpos_sur_IS(k+2*cand_num_max) - y_con2)**2
                        ind12 = k
                    else
                        if ( (x_IS - relpos_sur_IS(k) - x_con2)**2 + &
                             (y_IS - relpos_sur_IS(k+2*cand_num_max) - y_con2)**2 < dis12 ) then
                            ind12 = k
                        end if
                    end if
                end if
            end do
        else
            FP(1) = FP2
            pos_FP(1) = x2_reli
            pos_FP(1+FP_num_max) = y2_reli
            FP(2) = FP1
            pos_FP(2) = x1_reli
            pos_FP(2+FP_num_max) = y1_reli
            pos_FI(1) = x_IS
            pos_FI(1+FP_num_max) = y_IS
            start = FP2
            con1 = FP2
            con2 = FP1
            x_con2 = x1_reli
            y_con2 = y1_reli
            ind12 = ind_FP2(1)
        end if
        ! circle chain
        do j = 1, FP_num_max
            ! calc. next interseciton and particle
            ind_next = ind12 - 1
            if ( ind_next == 0 ) then
                ind_next = len_sur_part
            end if
            x_IS = relpos_sur_IS(ind_next) + x_con2
            y_IS = relpos_sur_IS(ind_next+2*cand_num_max) + y_con2
            pos_FI(j+1) = x_IS
            pos_FI(j+1+FP_num_max) = y_IS
            next_part = sur_part(ind_next)
            if ( next_part == start ) exit
            x_con2 = pos(next_part) - xi
            y_con2 = pos(next_part+num_part) - yi
            if ( x_con2 >= lx/2 ) then
                x_con2 = x_con2 - lx
            else if ( x_con2 < -lx/2 ) then
                x_con2 = x_con2 + lx
            end if
            if ( y_con2 >= ly/2 ) then
                y_con2 = y_con2 - ly
            else if ( y_con2 < -ly/2 ) then
                y_con2 = y_con2 + ly
            end if
            FP(j+2) = next_part
            pos_FP(j+2) = x_con2
            pos_FP(j+2+FP_num_max) = y_con2
            con1 = con2
            con2 = next_part
            ! calc. surrounding IS and the index of intersection
            call calc_sur_IS(con2,i,lx,ly,num_part,pos,rad_list,mask_num,gd_neigh_list,&
                             cand_num_max,radi,radi2,sur_part,relpos_sur_IS,len_sur_part)
            dis12 = -1.d0
            do k = 1, len_sur_part
                if ( sur_part(k) == con1 ) then
                    if ( dis12 == -1.d0 ) then
                        dis12 = (x_IS - relpos_sur_IS(k) - x_con2)**2 + &
                                (y_IS - relpos_sur_IS(k+2*cand_num_max) - y_con2)**2
                        ind12 = k
                    else
                        if ( (x_IS - relpos_sur_IS(k) - x_con2)**2 + &
                             (y_IS - relpos_sur_IS(k+2*cand_num_max) - y_con2)**2 < dis12 ) then
                            ind12 = k
                        end if
                        exit
                    end if
                end if
            end do
        end do
        num_FP = j + 1
        ! calc. FV, FS
        poly_area = 0.d0
        seg_area = 0.d0
        seg_len = 0.d0
        do j = 1, num_FP
            j_next = mod(j,num_FP) + 1
            poly_area = poly_area + pos_FI(j)*pos_FI(j_next+FP_num_max) &
                                  - pos_FI(j+FP_num_max)*pos_FI(j_next)
            x1_reli = pos_FI(j) - pos_FP(j_next)
            y1_reli = pos_FI(j+FP_num_max) - pos_FP(j_next+FP_num_max)
            x2_reli = pos_FI(j_next) - pos_FP(j_next)
            y2_reli = pos_FI(j_next+FP_num_max) - pos_FP(j_next+FP_num_max)
            angle = dabs(datan2(y1_reli,x1_reli) - datan2(y2_reli,x2_reli))
            if ( angle > Pi ) then
                angle = 2.d0 * Pi - angle
            end if
            seg_area = seg_area + angle * (rad_list(FP(j_next)) + radi)**2 &
                                - (-x1_reli * y2_reli + x2_reli * y1_reli)
            seg_len = seg_len + angle*(rad_list(FP(j_next)) + radi)
        end do
        FV(i) = (poly_area - seg_area) / 2.d0
        FS(i) = seg_len
    end do
    return
end subroutine

subroutine calc_sur_IS(target_part,exclude_part,lx,ly,num_part,pos,rad_list,mask_num,gd_neigh_list,&
                       cand_num_max,radi,radi2,sur_part,relpos_sur_IS,len_sur_part)
    implicit none
    integer, intent(in) :: target_part
    integer, intent(in) :: exclude_part
    double precision, intent(in) :: lx, ly
    integer, intent(in) :: num_part
    double precision, dimension(2*num_part), intent(in) :: pos
    double precision, dimension(num_part), intent(in) :: rad_list
    integer, intent(in) :: mask_num
    integer, dimension(mask_num,num_part), intent(in) :: gd_neigh_list
    integer, intent(in) :: cand_num_max
    double precision, intent(in) :: radi
    double precision, intent(in) :: radi2
    integer, dimension(2*cand_num_max), intent(out) :: sur_part
    double precision, dimension(4*cand_num_max), intent(out) :: relpos_sur_IS
    integer, intent(out) :: len_sur_part

    double precision :: dis_x, dis_y
    double precision :: dis
    double precision :: dis_j
    double precision, dimension(2*cand_num_max) :: ang_list
    double precision, dimension(4) :: IS ! InterSection
    integer j, k, m

    ! cutoff
    k=1
    do m = 1, gd_neigh_list(mask_num,target_part)
        j = gd_neigh_list(m,target_part)
        if ( j == exclude_part ) cycle
        if ( j /= 0 ) then
            dis_x = pos(j) - pos(target_part)
            dis_y = pos(j+num_part) - pos(target_part+num_part)
            if ( dis_x >= lx/2 ) then
                dis_x = dis_x - lx
            else if ( dis_x < -lx/2 ) then
                dis_x = dis_x + lx
            end if
            if ( dis_y >= ly/2 ) then
                dis_y = dis_y - ly
            else if ( dis_y < -ly/2 ) then
                dis_y = dis_y + ly
            end if 
            dis = dsqrt(dis_x**2 + dis_y**2)
            dis_j = dis - rad_list(target_part) - rad_list(j)
            if ( dis_j <= radi2 ) then
                sur_part(2*k-1) = j
                sur_part(2*k) = j
                call calc_IS(radi,dis,0.d0,0.d0,rad_list(target_part),dis_x,dis_y,rad_list(j),IS)
                relpos_sur_IS(2*k-1) = IS(1) 
                relpos_sur_IS(2*k-1+2*cand_num_max) = IS(2)
                relpos_sur_IS(2*k) = IS(3)
                relpos_sur_IS(2*k+2*cand_num_max) = IS(4)
                ang_list(2*k-1) = datan2(IS(2),IS(1))
                ang_list(2*k) = datan2(IS(4),IS(3))
                k = k + 1
            end if
        end if
    end do
    len_sur_part = 2 * (k - 1)
    call sort_integer_nochange(2*cand_num_max,len_sur_part,ang_list,sur_part)
    call sort_double_nochange(2*cand_num_max,len_sur_part,&
                              ang_list,relpos_sur_IS(:2*cand_num_max))
    call sort_double_nochange(2*cand_num_max,len_sur_part,&
                              ang_list,relpos_sur_IS(2*cand_num_max+1:4*cand_num_max))
    return
end subroutine

subroutine calc_IS(radi,dis12,pos_1x,pos_1y,rad_1,pos_2x,pos_2y,rad_2,IS)
    implicit none
    double precision, intent(in) :: radi
    double precision, intent(in) :: dis12
    double precision, intent(in) :: pos_1x, pos_1y
    double precision, intent(in) :: rad_1
    double precision, intent(in) :: pos_2x, pos_2y
    double precision, intent(in) :: rad_2
    double precision, dimension(4), intent(out) :: IS
    double precision :: center
    double precision :: dis_vert
    double precision, dimension(2) :: vec_cent, vec_vert1, vec_vert2

    center = ((rad_2 + radi)**2 - (rad_1 + radi)**2 + dis12**2)*0.5 / dis12
    vec_cent(1) = (pos_1x - pos_2x)/dis12 * center
    vec_cent(2) = (pos_1y - pos_2y)/dis12 * center
    dis_vert = dsqrt((rad_2 + radi)**2 - center**2)
    vec_vert1(1) = -vec_cent(2) / center
    vec_vert1(2) =  vec_cent(1) / center
    vec_vert2(1) =  vec_cent(2) / center
    vec_vert2(2) = -vec_cent(1) / center
    IS(1) = pos_2x + vec_cent(1) + dis_vert * vec_vert1(1)
    IS(2) = pos_2y + vec_cent(2) + dis_vert * vec_vert1(2)
    IS(3) = pos_2x + vec_cent(1) + dis_vert * vec_vert2(1)
    IS(4) = pos_2y + vec_cent(2) + dis_vert * vec_vert2(2)
    return
end subroutine
