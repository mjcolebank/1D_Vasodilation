subroutine impedance_driver(tmstps,Period,ff1,ff2,ff3,rho,mu,r_root,r_min,y_xt,Lr,Fr2,q,g,&
                            par1_in,par2_in,lrr,Z0,term_ID,Z_flag,vaso)
  use f90_tools
  use root_imp 
  implicit none

  integer, intent(in)    :: tmstps
  real(lng), intent(in)  :: Period,ff1,ff2,ff3,rho,mu,r_root,r_min,Lr,Fr2,q,g
  real(lng), intent(in)  :: par1_in,par2_in,lrr,Z0,vaso
  real(lng), intent(out) :: y_xt(tmstps)
  integer, intent(in)    :: term_ID, Z_flag


  call impedance (tmstps,Period,ff1,ff2,ff3,rho,mu,r_root,r_min,y_xt,Lr,Fr2,q,g,par1_in,par2_in,lrr,Z0, term_ID,Z_flag,vaso)
end subroutine impedance_driver
