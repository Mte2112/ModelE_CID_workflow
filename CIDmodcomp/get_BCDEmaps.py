def bcde_map(ds1, ds2):
    '''
    ds1 should be time-averaged map for model dataset 1
    ds2 should be time-averaged map for model dataset 2
    ds1 - ds2, i.e.
    B = model A map - default A map
    C = model A map - observations A map (make sure to interpolate model to obs grid. can use cdo remap)
    D = model A map - model pre-industrial A map
    E = model D map - default D map
    '''

    out_ds = ds1 - ds2

    return out_ds
