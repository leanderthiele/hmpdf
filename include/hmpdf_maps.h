#ifndef HMPDF_MAPS_H
#define HMPDF_MAPS_H

#include "hmpdf_object.h"

int hmpdf_get_map_op(hmpdf_obj *d,
                     int Nbins,
                     double binedges[Nbins+1],
                     double op[Nbins],
                     int new_map);

int hmpdf_get_map(hmpdf_obj *d,
                  double **map,
                  long *Nside,
                  int new_map);

#endif
