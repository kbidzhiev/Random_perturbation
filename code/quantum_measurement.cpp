#include "itensor/all.h"
#include "observables.hpp"

using namespace itensor;
using namespace std;








void FlipSpin(MPS& psi, const SiteSet &sites, const int j, Args args){
	auto l1 = Index(12,"l1");
	auto l2 = Index(12,"l2");
	auto s1 = Index(4,"s1");
	auto s2 = Index(4,"s2");

	auto H = randomITensorC(s1,s2,prime(s1),prime(s2));
	H = 0.5*(H+dag(swapTags(H,"0","1")));

	auto [Q,D] = diagHermitian(H);

	auto U = expHermitian(ITensor H);

	cout << U;

//	auto j = gate.i1;
//	auto &G = gate.G;
//	psi.position(j);
//	auto WF = psi(j) * psi(j + 1) * psi(j + 2);
//	WF = G * WF;
//	WF /= norm(WF);
//	WF.noPrime();
//	{
//		auto [Uj1, Vj1] = factor(WF,
//				{ siteIndex(psi, j), leftLinkIndex(psi, j) }, args);
//		auto indR = commonIndex(Uj1, Vj1);
//		auto [Uj2, Vj2] = factor(Vj1, { siteIndex(psi, j + 1), indR }, args);
//		psi.set(j, Uj1);
//		psi.set(j + 1, Uj2);
//		psi.set(j + 2, Vj2);
}
