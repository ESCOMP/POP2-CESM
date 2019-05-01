module POP_MCT_vars_mod

   use mct_mod
   use kinds_mod

   implicit none
   save
   public

   integer(int_kind) :: POP_MCT_OCNID
   type(mct_gsMap), pointer :: POP_MCT_gsMap_o     ! 2d, points to cdata
   type(mct_gGrid), pointer :: POP_MCT_dom_o       ! 2d, points to cdata
   type(mct_gsMap)          :: POP_MCT_gsMap3d_o   ! for 3d streams, local
   type(mct_gGrid)          :: POP_MCT_dom3d_o     ! for 3d streams, local

end module POP_MCT_vars_mod
