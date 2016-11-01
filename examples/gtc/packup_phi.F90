  subroutine packup_phi(ipack)
  use field_array
  use particle_array
  integer ipack

  if(ipack==1)then
    packphi(1:3,:,:)=gradphi(1:3,:,:)
    packphi(4,:,:)=0
  elseif(ipack==2)then
    packphi(4,:,:)=phit(:,:)
  else
    write(6,*)'ERROR:  Wrong type of packing of phi!'
    stop
  endif
end subroutine
  