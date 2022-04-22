subroutine plumes(V,L,AG,a,row,col,kend)
implicit none

!!January 2018 

! The routine identifies all four-connected elements in a matrix. 
! The returns are:
! V - arbitrary island number
! L - nx3 matrix where earch row has information about a particular cluster.
!     the first column is the cluster number which corresponds to V. The second
!     column is the number of elements in the cluster, and the third column is 
!     the value of the elements in that cluster (the elements of the input field).
! AG - Binary matrix with ones in the positions of the elements of the largers 
!      four-connected clusters. 

integer, intent(in) :: row,col,kend
integer, intent(in), dimension(row,col) :: a
integer, intent(out), dimension(kend) :: L,V
integer :: cnt,pp,mm,kk
integer :: i,j,cntr,gg,hh,ii,jj,idxx,IDX
integer, intent(out), dimension(row,col) :: AG



do j=1,col
 do i=1,row
  AG(i,j) = 0  !! This matrix will hold all of the islands found.
 enddo
enddo


do j=1,kend
V(j)=0
enddo

L = V  !! Hold the number of elements in each island.
cntr = 1  !! Label the individual islands by the order they are found.

!! Notes on comments in the loops:  CRNT is the element of the matrix we are
!! currently looking at.  RGHT is the element to the immediate right of
!! CRNT.  RUD is the element to the right and up from the CRNT.
!! The RUD element is directly above RGHT.

do gg = 1,col-1  !! Look along the first row.
    if (a(1,gg) == a(1,gg+1))then  !! CRNT matches RGHT.
        if (AG(1,gg) == 0)then  !! CRNT does not have an island number.
            AG(1,gg) = cntr  !! Assign an island number to CRNT.
            AG(1,gg+1) = cntr  !! Assign an island number to RGHT.
            V(cntr) = a(1,gg)  !! Assign the value of the island.
            L(cntr) = 2  !! Add to the island count.
            cntr = cntr + 1  !! Increment the counter.
        else  !! CRNT does have an island number.
            AG(1,gg+1) = AG(1,gg)  !! Assign the island number to RGHT.
            L(AG(1,gg)) = L(AG(1,gg)) + 1  !! Add to the island count.
        endif
    endif
enddo


do hh = 1,row-1  !! Look down the first column.
    if (a(hh,1)==a(hh+1,1))then  !! CRNT matches 'RGHT'.
        if (AG(hh,1)==0)then  !! CRNT does not have an island number.
            AG(hh,1) = cntr  !! Assign an island number to CRNT.
            AG(hh+1,1) = cntr  !! Assign an island number to 'RGHT'.
            V(cntr) = a(hh,1)  !! Assign the value of the island.
            L(cntr) = 2  !! Add to the island count.
            cntr = cntr + 1
        else  !! CRNT does have an island number.
            AG(hh+1,1) = AG(hh,1)  !! Assign an island number to 'RGHT'.
            L(AG(hh,1)) = L(AG(hh,1)) + 1  !! Add to the island count.
        endif
    endif
enddo


!! Now we can look at the rest of the matrix.
do ii = 2,row  !! Start on the second row.
    do jj = 1,col-1  !! Start on the first column.
        if (a(ii,jj)==a(ii,jj+1))then  !! CRNT matches RGHT.
            if (a(ii,jj)==a(ii-1,jj+1))then  !! CRNT matches RUD too.
                if  (AG(ii,jj)==0 .and. AG(ii-1,jj+1)==0)then  !! Both aren't yet grouped.
                    AG(ii,jj) = cntr  !! Give CRNT a new island number.
                    AG(ii-1,jj+1) = cntr  !! Give RUD the new island num.
                    AG(ii,jj+1) = cntr  !! Give RGHT the new island number.
                    V(cntr) = a(ii,jj) !! Store the value of the island.
                    L(cntr) = 3  !! Number of members in the new island.
                    cntr = cntr + 1  !! Increment the counter.
                elseif (AG(ii-1,jj+1)==0)then  !! RUD not yet been grouped, CRNT is.
                    AG(ii-1,jj+1) = AG(ii,jj) !! Give RUD CRNTs isl. num.
                    AG(ii,jj+1) = AG(ii,jj)  !! And RGHT as well.
                    L(AG(ii,jj)) = L(AG(ii,jj)) + 2  !! Add to island size.
                elseif (AG(ii,jj)==0)then  !! CRNT not yet grouped, RUD is grouped.
                    AG(ii,jj) = AG(ii-1,jj+1) !! Give CRNT RUDs isl. num.
                    AG(ii,jj+1) = AG(ii-1,jj+1)  !! And RGHT as well.
                    L(AG(ii-1,jj+1)) = L(AG(ii-1,jj+1)) + 2  !! Add to cnt.
                else !! Both CRNT and RUD have been grouped:  merge islands.
                    !! First decide which island has the least members.
                    if (L(AG(ii-1,jj+1))<L(AG(ii,jj)))then
                        idxx = AG(ii-1,jj+1)  !! Save RUDs island number.
                        IDX = AG(ii,jj) !! Used below: new islands.
                    else
                        idxx = AG(ii,jj)  !! Save CRNTs island number.
                        IDX = AG(ii-1,jj+1) !! Used below: new islands.
                    endif
                    cnt = 1 !! The counter.
                    if (idxx .ne. IDX)then !! They could already match!
                        !! These next do loops are faster than using find,
                        !! indexing the whole matrix.  We are searching do
                        !! the old island numbers to add to the updated isl.
                        do pp = ii+1,row !! Must search first column too.
                            if (AG(pp,1)==idxx)then
                                AG(pp,1) = IDX  !! Assign island number.
                                cnt = cnt + 1  !! Increment member counter.
                            else
                                goto 100 !! Stop search if one mismatch found.
                            endif
                        enddo
100                     do mm = ii,-1,1  !! Start at current row, work up.
                            do kk = 1,col
                                if (AG(mm,kk) == idxx)then
                                    AG(mm,kk) = IDX !! Assign new isl. num.
                                    cnt = cnt + 1 !! Increment count.
                                    if (cnt > L(idxx))then
                                        goto 101 !! Stop search do old islnds.
                                    endif
101                             endif
                            enddo
                            if (cnt > L(idxx))then
                                goto 102 !! Stop search do old island numbers.
102                         endif
                        enddo
                    endif
                    AG(ii,jj+1) = IDX  !! Give RGHT the island number.
                    L(IDX) = L(IDX) +  cnt  !! Add to count.
                    L(idxx) = L(idxx) - cnt + 1  !! subtract from old island.
                endif
            elseif (AG(ii,jj)==0)then  !! RGHT matches CRNT, not RUD. Need new isl.
                AG(ii,jj) = cntr  !! Give CRNT a new island number.
                AG(ii,jj+1) = cntr  !! Give RGHT the new island number.
                V(cntr) = a(ii,jj)  !! Store the new island number.
                L(cntr) = 2  !! The number of members in the new island.
                cntr = cntr + 1  !! Increment the counter.
            else !! RGHT matches CRNT, not RUD.  No new island needed.
                AG(ii,jj+1) = AG(ii,jj)  !! Give RGHT CRNTs island number.
                L(AG(ii,jj)) = L(AG(ii,jj)) + 1 !! Add to island count.
            endif
        elseif (a(ii,jj+1)==a(ii-1,jj+1))then !! RUD & RGHT match, not CRNT.
            if (AG(ii-1,jj+1)==0)then !! RUD has not yet been grouped.
                AG(ii,jj+1) = cntr  !! Give RGHT new island number.
                AG(ii-1,jj+1) = cntr !! Give RUD new island number.
                V(cntr) = a(ii,jj+1) !! Store the value of the island.
                L(cntr) = 2 !! Add to island count.
                cntr = cntr + 1
            else !! RUD is already part of a island.
                AG(ii,jj+1) = AG(ii-1,jj+1) !! Add RGHT to RUD's island.
                L(AG(ii-1,jj+1)) = L(AG(ii-1,jj+1)) + 1 !! Add island cnt.
            endif
        endif
    enddo
enddo


!    NSV = [(1:cntr-1)',L(1:cntr-1)',V(1:cntr-1)'] !!2nd output.
!    NSV(NSV(:,2)==0,:)=[] !! Clear empty islands.


end subroutine plumes
