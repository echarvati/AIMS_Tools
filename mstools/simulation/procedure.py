class Procedure:
    NPT = 'npt'
    NPT_MULTI = 'npt-multi'
    NVT_MULTI = 'nvt-multi'
    NVT_CV = 'nvt-cv'
    NVT_VISCOSITY = 'nvt-viscosity'
    NVT_VACUUM = 'nvt-vacuum'
    NVT_SLAB = 'nvt-slab'
    NPT_BINARY_SLAB = 'npt-binary-slab'
    NPT_PPM = 'ppm'
    choices = [NPT, NVT_CV, NVT_VISCOSITY, NVT_VACUUM, NVT_SLAB, NPT_BINARY_SLAB, NPT_PPM, NPT_MULTI, NVT_MULTI]

    prior = {
        NVT_CV : NPT,
        NVT_VISCOSITY: NPT,
        NPT_PPM: NPT,
        NVT_MULTI: NPT,
    }
