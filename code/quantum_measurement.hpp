#pragma once

#include "itensor/all.h"
#include <vector>
#include <string>
#include <complex>

using namespace itensor;
using namespace std;




void FlipSpin(MPS& psi, const SiteSet &sites, const int j, Args args);
