class Procedure:
    NPT = 'npt'
    NPT_MULTI = 'npt-multi'
    NPT_V_RESCALE = 'npt-v-rescale'
    NVT_MULTI = 'nvt-multi'
    NVT_CV = 'nvt-cv'
    NVT_VISCOSITY = 'nvt-viscosity'
    NVT_VACUUM = 'nvt-vacuum'
    NVT_SLAB = 'nvt-slab'
    NPT_BINARY_SLAB = 'npt-binary-slab'
    NPT_PPM = 'ppm'
    NPT_2 = 'npt-2'
    NPT_3 = 'npt-2'
    NVT_MULTI_2 = 'nvt-multi-2'
    NVT_MULTI_3 = 'nvt-multi-3'
    choices = [NPT, NVT_CV, NVT_VISCOSITY, NVT_VACUUM, NVT_SLAB, NPT_BINARY_SLAB, NPT_PPM, NPT_MULTI, NVT_MULTI,
               NPT_V_RESCALE, NPT_2, NPT_3, NVT_MULTI_2, NVT_MULTI_3]

    prior = {
        NVT_CV: NPT,
        NVT_VISCOSITY: NPT,
        NPT_PPM: NPT,
        NVT_MULTI: NPT,
        NPT_2: NPT,
        NPT_3: NPT,
        NVT_MULTI_2: NPT_2,
        NVT_MULTI_3: NPT_3,
    }
