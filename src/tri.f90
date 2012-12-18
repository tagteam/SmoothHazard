!========================          TRI         ====================
!
! subroutine qui prend un vecteur, sa taille , et renvoit un vecteur trie et 
!       sa moyenne
       
       
       subroutine tri(vecteur,taille,moyenne)
       
       implicit none
       
       integer,intent(in)::taille
       double precision,dimension(taille),intent(inout)::vecteur
       double precision,intent(out)::moyenne
       double precision::temp
       integer::i,j
              
          i=1
          do while (i==1)
	     i=0
	     do j=1,taille-1
                if (vecteur(j) > vecteur(j+1)) then
	           temp=vecteur(j)
	           vecteur(j)=vecteur(j+1)
	           vecteur(j+1)=temp
	           i=1
	         end if  
	     end do     
          end do
	  
	  moyenne=0.d0
	  do i=1,taille
             moyenne=moyenne+vecteur(i)    
	  end do 
	  moyenne=moyenne/taille
       end subroutine tri
