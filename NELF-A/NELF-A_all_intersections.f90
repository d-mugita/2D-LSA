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
    double precision :: rad_M
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

    rad_M = (rad_S + rad_L) / 2.d0

!================ Main
    call gridmapping(lx,ly,num_part,pos,gd_num_x,gd_num_y,gd_map,gd_coord)
    call finding_neighbors(num_part,gd_num_x,gd_num_y,gd_map,gd_coord,mask_num,mask,gd_neigh_list)

    call NELFA_all_IS(Pi,lx,ly,num_part,rad_S,rad_L,pos,rad_list,mask_num,gd_neigh_list,&
                      cand_num_max,FP_num_max,rad_M,FV,FS)

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

subroutine NELFA_all_IS(Pi,lx,ly,num_part,rad_S,rad_L,pos,rad_list,mask_num,gd_neigh_list,&
                        cand_num_max,FP_num_max,rad_M,FV,FS)
    implicit none
    double precision, intent(in) :: Pi
    double precision, intent(in) :: lx, ly
    integer, intent(in) :: num_part
    double precision, intent(in) :: rad_S, rad_L
    double precision, dimension(2*num_part), intent(in) :: pos
    double precision, dimension(num_part), intent(in) :: rad_list
    integer, intent(in) :: mask_num
    integer, dimension(mask_num,num_part), intent(in) :: gd_neigh_list
    integer, intent(in) :: cand_num_max
    integer, intent(in) :: FP_num_max
    double precision, intent(in) :: rad_M
    double precision, dimension(num_part), intent(out) :: FV, FS

    ! AI, AP: All Intersections and their pair Particles
    ! FP, FI: Free volume construction Particles and Intersections
    double precision :: rad_S2, rad_L2
    double precision :: dis_x, dis_y
    double precision :: dis, dis_actual
    integer, dimension(cand_num_max) :: pair_list_S, pair_list_L
    double precision, dimension(2*cand_num_max) :: relpos_pair_S, relpos_pair_L
    double precision, dimension(cand_num_max) :: dis_list_S, dis_list_L
    integer, dimension(1) :: ind_min
    integer, dimension(num_part) :: FP1_list
    double precision, dimension(2*num_part) :: pos_FP1
    double precision, dimension(4) :: IS ! InterSection
    integer, dimension(num_part,2*cand_num_max) :: AP_S, AP_L
    double precision, dimension(num_part,4*cand_num_max) :: pos_AI_S, pos_AI_L
    double precision, dimension(2*cand_num_max) :: angle_list_S, angle_list_L
    double precision :: xi, yi
    double precision :: radi
    integer, dimension(num_part,2*cand_num_max) :: AP
    double precision, dimension(num_part,4*cand_num_max) :: pos_AI
    integer :: FP1, FP2
    double precision :: x1_reli, y1_reli, x2_reli, y2_reli
    integer :: len
    double precision, dimension(2*cand_num_max) :: dis2_i2IS
    integer, dimension(1) :: ind_FP2
    double precision :: x_IS, y_IS
    double precision :: vec1x, vec1y, vec2x, vec2y
    integer, dimension(FP_num_max) :: FP
    double precision, dimension(2*FP_num_max) :: pos_FP
    double precision, dimension(2*FP_num_max) :: pos_FI
    integer :: start
    integer :: con1, con2
    double precision :: x_con2, y_con2
    double precision :: dis12
    integer :: ind12, ind_next
    integer :: next_part
    integer :: num_FP
    double precision :: poly_area, seg_area, seg_len
    integer :: j_next
    double precision :: angle
    integer i, j, k, l, m

    ! cut-off 
    rad_S2 = 2.d0 * rad_S
    rad_L2 = 2.d0 * rad_L
    do i = 1, num_part
        k = 1
        l = 1
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
                dis = dsqrt(dis_x**2 + dis_y**2)
                dis_actual = dis - rad_list(i) - rad_list(j)
                if ( dis_actual <= rad_L2 ) then
                    pair_list_L(l) = j
                    relpos_pair_L(l) = dis_x
                    relpos_pair_L(l+cand_num_max) = dis_y
                    dis_list_L(l) = dis
                    l = l + 1
                    if ( dis_actual <= rad_S2 ) then
                        pair_list_S(k) = j
                        relpos_pair_S(k) = dis_x
                        relpos_pair_S(k+cand_num_max) = dis_y
                        dis_list_S(k) = dis - rad_list(j)
                        k = k + 1
                    end if
                end if
            end if
        end do
        k = k - 1
        l = l - 1
        ! calc. FP1
        ind_min = minloc(dis_list_S(:k))
        FP1_list(i) = pair_list_S(ind_min(1))
        pos_FP1(i) = relpos_pair_S(ind_min(1))
        pos_FP1(i+num_part) = relpos_pair_S(ind_min(1)+cand_num_max)
        ! calc. and sort AP
        do m = 1, k
            j = pair_list_S(m)
            AP_S(i,2*m-1) = j
            AP_S(i,2*m) = j
            call calc_IS(rad_S,dis_list_S(m)+rad_list(j),0.d0,0.d0,rad_list(i),&
                         relpos_pair_S(m),relpos_pair_S(m+cand_num_max),rad_list(j),IS)
            pos_AI_S(i,2*m-1) = IS(1) 
            pos_AI_S(i,2*m-1+2*cand_num_max) = IS(2)
            pos_AI_S(i,2*m) = IS(3)
            pos_AI_S(i,2*m+2*cand_num_max) = IS(4)
            angle_list_S(2*m-1) = datan2(IS(2),IS(1))
            angle_list_S(2*m) = datan2(IS(4),IS(3))
        end do
        AP_S(i,2*cand_num_max) = 2 * k
        call sort_integer_nochange(2*cand_num_max,2*k,angle_list_S,AP_S(i,:))
        call sort_double_nochange(2*cand_num_max,2*k,angle_list_S,pos_AI_S(i,:2*cand_num_max))
        call sort_double_nochange(2*cand_num_max,2*k,&
                                  angle_list_S,pos_AI_S(i,2*cand_num_max+1:4*cand_num_max))
        do m = 1 , l
            j = pair_list_L(m)
            AP_L(i,2*m-1) = j
            AP_L(i,2*m) = j
            call calc_IS(rad_L,dis_list_L(m),0.d0,0.d0,rad_list(i),&
                         relpos_pair_L(m),relpos_pair_L(m+cand_num_max),rad_list(j),IS)
            pos_AI_L(i,2*m-1) = IS(1) 
            pos_AI_L(i,2*m-1+2*cand_num_max) = IS(2)
            pos_AI_L(i,2*m) = IS(3)
            pos_AI_L(i,2*m+2*cand_num_max) = IS(4)
            angle_list_L(2*m-1) = datan2(IS(2),IS(1))
            angle_list_L(2*m) = datan2(IS(4),IS(3))
        end do
        AP_L(i,2*cand_num_max) = 2 * l
        call sort_integer_nochange(2*cand_num_max,2*l,angle_list_L,AP_L(i,:))
        call sort_double_nochange(2*cand_num_max,2*l,angle_list_L,pos_AI_L(i,:2*cand_num_max))
        call sort_double_nochange(2*cand_num_max,2*l,&
                                  angle_list_L,pos_AI_L(i,2*cand_num_max+1:4*cand_num_max))
    end do
    ! calc. pos_FP and pos_FI
    do i = 1, num_part
        xi = pos(i)
        yi = pos(i+num_part)
        radi = rad_list(i)
        if ( radi < rad_M ) then
            AP = AP_S
            pos_AI = pos_AI_S
        else
            AP = AP_L
            pos_AI = pos_AI_L
        end if
        ! FP1
        FP1 = FP1_list(i)
        x1_reli = pos_FP1(i)
        y1_reli = pos_FP1(i+num_part)
        ! FP2
        len = AP(FP1,2*cand_num_max)
        do j = 1, len
            dis2_i2IS(j) = (pos_AI(FP1,j) + x1_reli)**2 + (pos_AI(FP1,j+2*cand_num_max) + y1_reli)**2
        end do
        ind_FP2 = minloc(dis2_i2IS(:len))
        if ( AP(FP1,ind_FP2(1)) == i ) then
            ind_FP2 = minloc(dis2_i2IS(:len), mask = dis2_i2IS(:len) > dis2_i2IS(ind_FP2(1)))
            if ( AP(FP1,ind_FP2(1)) == i ) then
                ind_FP2 = minloc(dis2_i2IS(:len), mask = dis2_i2IS(:len) > dis2_i2IS(ind_FP2(1)))
            end if
        end if
        FP2 = AP(FP1,ind_FP2(1))
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
        x_IS = pos_AI(FP1,ind_FP2(1)) + x1_reli
        y_IS = pos_AI(FP1,ind_FP2(1)+2*cand_num_max) + y1_reli
        ! decide con1, con2 (con: confirmed particle)
        vec1x = x1_reli - x_IS ! vec1: vector from IS to FP1
        vec1y = y1_reli - y_IS
        vec2x = x2_reli - x_IS ! vec2: vector from IS to FP2
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
        end if
        ! circle chain
        do j = 1, FP_num_max
            ! calc. the index of intersection
            dis12 = -1.d0
            do k = 1, AP(con2,2*cand_num_max)
                if ( AP(con2,k) == con1 ) then
                    if ( dis12 == -1.d0 ) then
                        dis12 = (x_IS - pos_AI(con2,k) - x_con2)**2 + &
                                (y_IS - pos_AI(con2,k+2*cand_num_max) - y_con2)**2
                        ind12 = k
                    else
                        if( (x_IS - pos_AI(con2,k) - x_con2)**2 + &
                            (y_IS - pos_AI(con2,k+2*cand_num_max) - y_con2)**2 < dis12 ) then
                            ind12 = k
                        end if
                        exit
                    end if
                end if
            end do
            ! calc. next interseciton and particle
            ind_next = ind12 - 1
            if ( ind_next == 0 ) then
                ind_next = AP(con2,2*cand_num_max)
            end if
            if ( AP(con2,ind_next) == i ) then
                ind_next = ind_next - 1
                if ( ind_next == 0 ) then
                    ind_next = AP(con2,2*cand_num_max)
                end if
                if ( AP(con2,ind_next) == i ) then
                    ind_next = ind_next - 1
                    if ( ind_next == 0 ) then
                        ind_next = AP(con2,2*cand_num_max)
                    end if
                end if
            end if
            x_IS = pos_AI(con2,ind_next) + x_con2
            y_IS = pos_AI(con2,ind_next+2*cand_num_max) + y_con2
            pos_FI(j+1) = x_IS
            pos_FI(j+1+FP_num_max) = y_IS
            next_part = AP(con2,ind_next)
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
        end do
        num_FP = j + 1
        ! calc. FV, FS
        poly_area = 0.d0
        seg_area = 0.d0
        seg_len = 0.d0
        do j = 1, num_FP
            j_next = mod(j,num_FP) + 1
            poly_area = poly_area + pos_FI(j) * pos_FI(j_next+FP_num_max) &
                                  - pos_FI(j+FP_num_max) * pos_FI(j_next)
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
            seg_len = seg_len + angle * (rad_list(FP(j_next)) + radi)
        end do
        FV(i) = (poly_area - seg_area) / 2.d0
        FS(i) = seg_len
    end do
    return
end

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
