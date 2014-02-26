
       if(region_CPML(ispec)==1)then
! element is in the left cpml layer
                                              :: should all become CPML_X_ONLY = 1 in the new numbering
       else if(region_CPML(ispec)==2)then
! element is in the right cpml layer

======================================================================

       else if(region_CPML(ispec)==4)then
! element is in the top cpml layer
                                              :: should all become CPML_Z_ONLY = 2 in the new numbering
       else if(region_CPML(ispec)==8)then
! element is in the bottom cpml layer

======================================================================

       else if(region_CPML(ispec)==5)then
! element is in the left-top cpml corner

       else if(region_CPML(ispec)==6)then
! element is in the right-top cpml corner
                                              :: should all become CPML_XZ_ONLY = 3 in the new numbering
       else if(region_CPML(ispec)==9)then
! element is in the left-bottom cpml corner

       else if(region_CPML(ispec)==10)then
! element is in the right-bottom cpml corner

