#include "itensor/all.h"
#include "observables.hpp"

using namespace itensor;
using namespace std;






void FlipSpin(MPS& psi, const SiteSet &sites, const int j, Args args){
	auto j0 = siteIndex(psi, j);
	auto j1 = siteIndex(psi, j+1);
	auto T = randomITensor(j0,j1,prime(j0),prime(j1));
	auto W = T + swapTags(dag(T),"0","1");	//Make Hermitian matrix out of T
	auto U = expHermitian(W,1_i);	// compute exp(i * W) //	PrintData(U);
	psi.position(j);
	auto WF = psi(j) * psi(j + 1);
	WF = U * WF;
	WF /= norm(WF);
	WF.noPrime();
	auto [Uj1, Vj1] = factor(WF, { siteIndex(psi, j), leftLinkIndex(psi, j) }, args);
	psi.set(j, Uj1);
	psi.set(j + 1, Vj1);
}
