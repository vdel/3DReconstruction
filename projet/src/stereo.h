#ifndef ENERGYH
#define ENERGYH

#include "main.h"
#include "camera.h"
#include "carac.h"

double Ep(camvec& cl, const vec& p_world, Carac& CaracEstimate, double OutRatio = 0.5, double MinPb = 0.99);
double En(const SearchSpace& ss, Image<vector<Carac> >& carac, Image<pair<int, int> > depth_limit, const vec step, Coords<3> M);

void graphcut_optimize(bool verbose, SearchSpace& ss, int D, Image<pair<int, int> >& depth_limit, Image<int>& depth, Image<Carac>& caracs);
void lsqr_optimize(bool verbose, SearchSpace& ss, int D, Image<int>& depth, Image<Carac>& caracs, Image<double>& optimal_depth);

void show_energy_depth(SearchSpace& ss, Image<vector<double> > E, Image<pair<int, int> >& depth_limit, int D);
#endif
