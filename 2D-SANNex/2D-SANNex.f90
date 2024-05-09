program SANNex
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
    integer, parameter                     :: mask_num = 292 ! grid length 5: 120, 9: 292
    integer, dimension(2,mask_num)         :: mask
    integer, dimension(mask_num,num_part)  :: gd_neigh_list

!---- SANN
    ! effective diameter: The average distance between particles per collision
    double precision, parameter :: eff_dia = 0.964481663872923d0 ! 0.720
    !double precision, parameter :: eff_dia = 0.990731081374231d0 ! 0.760
    !double precision, parameter :: eff_dia = 1.003058806502530d0 ! 0.780
    ! non-dimensional cutoff radius: First minimun of radial distribution function)
    double precision, parameter :: cutoff_rad_nonD = 1.49303d0 ! 0.720
    !double precision, parameter :: cutoff_rad_nonD = 1.45347d0 ! 0.760 
    !double precision, parameter :: cutoff_rad_nonD = 1.43561d0 ! 0.780
    double precision            :: cutoff_rad
    double precision            :: cutoff_rad_cand
    integer, parameter          :: neigh_num_max = 120 ! bin, 3rd
    integer(8), dimension(num_part,neigh_num_max) :: neigh1
    integer, dimension(num_part)               :: neigh1_num
    integer(8), dimension(num_part,neigh_num_max) :: neigh2
    integer, dimension(num_part)               :: neigh2_num
    integer(8), dimension(num_part,neigh_num_max) :: neigh3
    integer, dimension(num_part)               :: neigh3_num

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
    call gridmask9(mask)

    cutoff_rad = cutoff_rad_nonD * eff_dia
    cutoff_rad_cand = 3.d0 * cutoff_rad ! grid length 5: 2.d0 *, 9: 3.d0 *

!================ Main
    call gridmapping(lx,ly,num_part,pos,gd_num_x,gd_num_y,gd_map,gd_coord)
    call finding_neighbors(num_part,gd_num_x,gd_num_y,gd_map,gd_coord,mask_num,mask,gd_neigh_list)

    neigh1 = 0
    neigh1_num = 0
    neigh2 = 0
    neigh2_num = 0
    neigh3 = 0
    neigh3_num = 0
    call SANN(Pi,lx,ly,num_part,pos,mask_num,gd_neigh_list,cutoff_rad_cand,neigh_num_max,&
              neigh1,neigh1_num,neigh2,neigh2_num,neigh3,neigh3_num)

    open(20, file='./neigh1.dat')
    open(21, file='./neigh1_num.dat')
    open(22, file='./neigh2.dat')
    open(23, file='./neigh2_num.dat')
    open(24, file='./neigh3.dat')
    open(25, file='./neigh3_num.dat')
    do i = 1, num_part
        write(20,*) neigh1(i,:) 
        write(21,*) neigh1_num(i)
        write(22,*) neigh2(i,:) 
        write(23,*) neigh2_num(i)
        write(24,*) neigh3(i,:) 
        write(25,*) neigh3_num(i)
    end do
    close(20)
    close(21)
    close(22)
    close(23)
    close(24)
    close(25)
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

subroutine sort(N,M,A,B)
    integer, intent(in) :: N ! array size
    integer, intent(in) :: M ! sorting part
    double precision, dimension(N), intent(inout) :: A ! array
    integer(8), dimension(N), intent(inout) :: B ! array
    double precision :: tmp1
    integer(8) :: tmp2
    integer :: i, j, k

    do i = 2, M
        tmp1 = A(i)
        tmp2 = B(i)
        do j = i-1, 1, -1
            if ( A(j) > tmp1 ) then
                A(j+1) = A(j)
                B(j+1) = B(j)
                k = j
            else
                k = j + 1
                exit
            end if
        end do
        A(k) = tmp1
        B(k) = tmp2
    end do
    return
end subroutine
    
subroutine SANN(Pi,lx,ly,num_part,pos,mask_num,gd_neigh_list,cutoff_rad_cand,neigh_num_max,&
                neigh1,neigh1_num,neigh2,neigh2_num,neigh3,neigh3_num)
    implicit none
    double precision, intent(in) :: Pi
    double precision, intent(in) :: lx, ly
    integer, intent(in) :: num_part
    double precision, dimension(2*num_part), intent(in) :: pos
    integer, intent(in) :: mask_num
    integer, dimension(mask_num,num_part), intent(in) :: gd_neigh_list
    double precision, intent(in) :: cutoff_rad_cand
    integer, intent(in) :: neigh_num_max
    integer(8), dimension(num_part,neigh_num_max), intent(out) :: neigh1
    integer, dimension(num_part), intent(out) :: neigh1_num
    integer(8), dimension(num_part,neigh_num_max), intent(out) :: neigh2
    integer, dimension(num_part), intent(out) :: neigh2_num
    integer(8), dimension(num_part,neigh_num_max), intent(out) :: neigh3
    integer, dimension(num_part), intent(out) :: neigh3_num

    integer :: num_cand
    integer(8), dimension(neigh_num_max) :: cand_NN
    double precision, dimension(neigh_num_max) :: dis_list  
    double precision :: dis_x, dis_y
    double precision :: dis2
    double precision :: sum_acos
    double precision :: cutoff_rad_SANN
    integer :: i, j, k, l, m, n    
    
    do i = 1, num_part
        dis_list = 0.d0
        k = 1
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
                dis2 = dis_x**2 + dis_y**2
                if ( dis2 <= cutoff_rad_cand**2 ) then
                    cand_NN(k) = j
                    dis_list(k) = dsqrt(dis2)
                    k = k + 1
                end if
            end if
        end do
        num_cand = k - 1
        call sort(neigh_num_max,num_cand,dis_list,cand_NN)

        ! first
        do l = 3, num_cand
            sum_acos = 0.d0
            cutoff_rad_SANN = dis_list(l+1)
            do j = 1, l
                sum_acos = sum_acos + dacos(dis_list(j)/cutoff_rad_SANN)
            end do
            if ( sum_acos >= Pi ) exit
        end do
        neigh1_num(i) = l
        do k = 1, l
            neigh1(i,k) = cand_NN(k)
        end do
        ! second
        do m = l+3, num_cand
            sum_acos = 0.d0
            cutoff_rad_SANN = dis_list(m+1)
            do j = l+1, m
                sum_acos = sum_acos + dacos(dis_list(j)/cutoff_rad_SANN)
            end do
            if ( sum_acos >= 2.d0*Pi ) exit
        end do
        neigh2_num(i) = m - l
        do k = l+1, m
            neigh2(i,k-l) = cand_NN(k)
        end do
        ! third
        do n = m+3, num_cand
            sum_acos = 0.d0
            cutoff_rad_SANN = dis_list(n+1)
            do j = m+1, n
                sum_acos = sum_acos + dacos(dis_list(j)/cutoff_rad_SANN)
            end do
            if ( sum_acos >= 3.d0*Pi ) exit
        end do
        neigh3_num(i) = n - m
        do k = m+1, n
            neigh3(i,k-m) = cand_NN(k)
        end do
    end do
    return
end subroutine
