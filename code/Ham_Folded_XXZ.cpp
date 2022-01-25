#include "itensor/all.h"
#include "time_evolution.hpp"
#include <cmath>

using namespace itensor;
using namespace std;


//
ThreeSiteHamiltonian::ThreeSiteHamiltonian(const SiteSet &sites,
		const ThreeSiteParam &param)
			: ampo(sites)
			, N(length(sites)) {
	//N = length(sites); // size of Hamiltonian comes with onject SiteSet sites, not just as number N
	init(param);   // initializing the Hamiltonian
	cout << "A Hamiltonian with " << N << " sites was constructed." << endl;
}

//Initialize Hamiltonian parameters
void ThreeSiteHamiltonian::init(const ThreeSiteParam &param) {
	const double J = param.val("J");
	double mu = 0;
	const double hL = param.val("hL");
	const double hR = param.val("hR");
	dot = N / 2 + 1;  //Position of the "dot"
	cout << "The dot is on site #" << dot << endl;
	for (int j = 1; j < N -1; ++j) {
		//Strange coefficients are needed to match
		// spin matrices and Pauli matrices -> Pauli = 2*Spin, so each matrix gives factor 2
		// one of 1/2 comes from the projector (1-sigma_z)/2
		// and the other is "Jacobian", i.e. 0.5 (SpSm+ SmSp) = SxSx + SySy
		ampo += J *  4 * 0.25, "S+", j, "S-", j + 2; //
		ampo += J *  4 * 0.25, "S-", j, "S+", j + 2;
		ampo += J * -8 * 0.25, "S+", j, "Sz", j + 1, "S-", j + 2;
		ampo += J * -8 * 0.25, "S-", j, "Sz", j + 1, "S+", j + 2;


		if (j <= dot) {
			mu = hL;
		} else {
			mu = hR;
		}
		ampo += mu, "Sz", j;
	}
	// boundary terms. for loop doesn't reach j == N-1 OR N
	ampo += mu, "Sz", N - 1;
	ampo += mu, "Sz", N;
}




// Folded XYZ with integrability breaking term
//https://arxiv.org/abs/2110.11322 in SIGMA repr Eq(1)
HamiltonianFoldedXYZ::HamiltonianFoldedXYZ(const SiteSet &sites,
		const ThreeSiteParam &param)
			: ampo(sites)
			, N(length(sites)) {
	//N = length(sites); // size of Hamiltonian comes with onject SiteSet sites, not just as number N
	init(param);   // initializing the Hamiltonian
	cout << "A Hamiltonian with " << N << " sites was constructed." << endl;
}

//Initialize Hamiltonian parameters
void HamiltonianFoldedXYZ::init(const ThreeSiteParam &param) {
	const double Jx = param.val("J");
	const double Jy = param.val("Jy");
	const double Delta = param.val("Delta");
	const double J2 = param.val("J2");

	for (int j = 1; j < N - 1; ++j) {
		//XYZ in sigma basis
		ampo += Jx * 4, "Sx", j, "Sx", j + 2;
		ampo += -Jy * 8, "Sx", j, "Sz", j + 1, "Sx", j + 2;
		ampo += Delta * 2, "Sz", j;
		// integrability breaking term
		ampo += J2 * 4, "Sz", j, "Sz", j + 1;
	}
	// boundary terms. "for" loop doesn't reach j == N-1 and N
	ampo += Delta * 2, "Sz", N - 1;
	ampo += Delta * 2, "Sz", N;
	ampo += J2 * 4, "Sz", N-1, "Sz", N;
}


