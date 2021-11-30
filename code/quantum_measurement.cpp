#include "itensor/all.h"
#include "observables.hpp"

using namespace itensor;
using namespace std;








void FlipSpin(MPS& psi, const SiteSet &sites, const int j, Args args){
//	auto l1 = Index(12,"l1");
//	auto l2 = Index(12,"l2");
//	auto s1 = Index(4,"s1");
//	auto s2 = Index(4,"s2");
//
//	auto H = randomITensorC(s1,s2,prime(s1),prime(s2));
//	H = 0.5*(H+dag(swapTags(H,"0","1")));
//
//	auto [Q,D] = diagHermitian(H);
//
//	auto U = expHermitian(H);
//
//	cout << U << endl;
//	Print(D);


	auto j0 = siteIndex(psi, j);
	auto j1 = siteIndex(psi, j+1);

	auto T = randomITensor(j0,j1,prime(j0),prime(j1));

	//Make Hermitian matrix out of T
	auto W = T + swapTags(dag(T),"0","1");


	auto U = expHermitian(W,1_i);	// compute exp(i * W)
//	cout << U << endl;
//	PrintData(U);
	psi.position(j);
	auto WF = psi(j) * psi(j + 1);
	WF = U * WF;
	WF /= norm(WF);
	WF.noPrime();
	auto [Uj1, Vj1] = factor(WF, { siteIndex(psi, j), leftLinkIndex(psi, j) }, args);
	psi.set(j, Uj1);
	psi.set(j + 1, Vj1);

}
