#include "itensor/all.h"
#include "time_evolution.hpp"
#include <cmath>

using namespace itensor;
using namespace std;


//Trotter Gates
TrotterExp::TrotterExp(const SiteSet &sites, const ThreeSiteParam &param,
		const complex<double> tau) {
	initialize(sites, param, tau);
}
void TrotterExp::initialize(const SiteSet &sites, const ThreeSiteParam &param,
		const complex<double> tau) {
	//const int begin = param.val("begin");
	const int begin = 1;
	const int end = param.val("N");
	const int order = param.val("TrotterOrder");
	if (order == 1) {
		cout << "trotter 1 scheme" << endl;

		TimeGates(begin, end, tau, sites, param);
		TimeGates(begin + 1, end, tau, sites, param);
		TimeGates(begin + 2, end, tau, sites, param);

	} else {
		cout << "trotter 2 scheme" << endl;
		/*
		 double a1 = 1. / 6;		// more precise arrpoximation coefficients
		 double a2 = 1 - 2. * a1;
		 double b1 = (3 - sqrt(3)) / 6.;
		 double b2 = 1. / 2 - b1;
		 double c1 = 1. / 2;
		 */
		double begin0 = begin; //this variable are needed to change operators ABC
		double begin2 = begin + 1;
		double begin4 = begin + 2;
		//Trotter gates from arxiv.org/abs/1901.04974
		// Eq. (38),(47)

		cout << "Time evolutions " << endl;
		TimeGates(begin0, end, 0.5 * tau, sites, param); //A
		TimeGates(begin2, end, 0.5 * tau, sites, param); //B
		TimeGates(begin4, end, tau, sites, param); //C
		TimeGates(begin2, end, 0.5 * tau, sites, param); //B
		TimeGates(begin0, end, 0.5 * tau, sites, param); //A
		/*
		 TimeGates(begin0, end, a1 * tau, sites, param); //A
		 TimeGates(begin2, end, b1 * tau, sites, param); //B
		 TimeGates(begin4, end, c1 * tau, sites, param); //C
		 TimeGates(begin2, end, b2 * tau, sites, param); //B
		 TimeGates(begin0, end, a2 * tau, sites, param); //A
		 TimeGates(begin2, end, b2 * tau, sites, param); //B
		 TimeGates(begin4, end, c1 * tau, sites, param); //C
		 TimeGates(begin2, end, b1 * tau, sites, param); //B
		 TimeGates(begin0, end, a1 * tau, sites, param); //A
		 */
	}
}

void TrotterExp::TimeGates(const int begin, const int end,
		const complex<double> tau, const SiteSet &sites,
		const ThreeSiteParam &param) {
	const int step = 3;
	const double J = param.val("J");
	//cout << "Gates starts from " << begin << endl;
	for (int j = begin; j < end - 1; j += step) {
		//cout << "j = (" << j << ", " << j + 1 << ", " << j + 2 << ")"
		//		<< endl;
		//this part act on real sites
		auto hh = J * 4 * 0.25 * op(sites, "Sp", j) * op(sites, "Id", j + 1)
				* op(sites, "Sm", j + 2);
		hh += J * 4 * 0.25 * op(sites, "Sm", j) * op(sites, "Id", j + 1)
				* op(sites, "Sp", j + 2);
		hh += -J * 8 * 0.25 * op(sites, "Sp", j) * op(sites, "Sz", j + 1)
				* op(sites, "Sm", j + 2);
		hh += -J * 8 * 0.25 * op(sites, "Sm", j) * op(sites, "Sz", j + 1)
				* op(sites, "Sp", j + 2);



		auto G = expHermitian(hh, tau);
		gates.emplace_back(j, move(G));
	}
}
void TrotterExp::Evolve(MPS &psi, const Args &args) {
	for (auto &gate : gates) {
		auto j = gate.i1;
		auto &G = gate.G;
		psi.position(j);
		auto WF = psi(j) * psi(j + 1) * psi(j + 2);
		WF = G * WF;
		WF /= norm(WF);
		WF.noPrime();
		{
			auto [Uj1, Vj1] = factor(WF,
					{ siteIndex(psi, j), leftLinkIndex(psi, j) }, args);
			auto indR = commonIndex(Uj1, Vj1);
			auto [Uj2, Vj2] = factor(Vj1, { siteIndex(psi, j + 1), indR }, args);
			psi.set(j, Uj1); // instead of copying I can move here??
			psi.set(j + 1, Uj2);
			psi.set(j + 2, Vj2);

		}
	}
}

vector<MPO> XXZ_time_evol(const SiteSet &sites, const ThreeSiteParam &param) {
	double tau = param.val("tau");
	const int order = param.val("TrotterOrderXXZ");
	vector<complex<double>> time_steps;
	time_steps.reserve(7);
	if (order == 1) {
		//Approx. with error O(tau^3)
		time_steps.push_back(Cplx_i * tau);
	}
	if (order == 2) {
		//Approx. with error O(tau^3)
		time_steps.push_back(0.5 * ( 1 + Cplx_i) * tau);
		time_steps.push_back(0.5 * (-1 + Cplx_i) * tau);
	}
	if (order == 3) {
		//Approx. with error O(tau^4) [thanks to Kemal]
		//Tested -it's ok
		time_steps.push_back(0.10566243270259355887 - 0.39433756729740644113 * Cplx_i);
		time_steps.push_back(Cplx_i * time_steps[0]);
		time_steps.push_back(conj(time_steps[1]));
		time_steps.push_back(Cplx_i * time_steps[2]);
		for (auto &time_step : time_steps){
			time_step *= Cplx_i * tau;
		}
	}
	if (order == 4) {
		//Approx. with error O(tau^5) [thanks to Kemal twice]
		//WARNINIG Not fully TESTED
		//Tested -it's ok
		time_steps.push_back( 0.25885339861091821723 + 0.04475613401114190287 * Cplx_i);
		time_steps.push_back(-0.03154685814880379274 + 0.24911905427556321757 * Cplx_i);
		time_steps.push_back( 0.19082905211066719664 - 0.23185374923210605447 * Cplx_i);
		time_steps.push_back( 0.1637288148544367438753);
		time_steps.push_back(conj(time_steps[2]));
		time_steps.push_back(conj(time_steps[1]));
		time_steps.push_back(conj(time_steps[0]));
		for (auto &time_step : time_steps){
			time_step *= Cplx_i * tau;
		}

	}

	vector<MPO> Exp_H_vec;
	Exp_H_vec.reserve(time_steps.size());

	XXZ H_XXZ(sites, param);
	for (auto time_step : time_steps){
		MPO expH_i = toExpH(H_XXZ.ampo, time_step);
		Exp_H_vec.push_back(expH_i);
	}

	return Exp_H_vec;
}



///////////////////////////////////////////////////////
//XXZ trotter gates

//Trotter Gates
TrotterExpXXZ::TrotterExpXXZ(const SiteSet &sites, const ThreeSiteParam &param,
		const complex<double> tau) {
	initialize(sites, param, tau);
}
void TrotterExpXXZ::initialize(const SiteSet &sites, const ThreeSiteParam &param,
		const complex<double> tau) {
	//const int begin = param.val("begin");
	const int begin = 1;
	const int end = param.val("N");
	const int order = param.val("TrotterOrder");
	if (order == 1) {
		cout << "trotter 1 scheme" << endl;
		TimeGates(begin, end, tau, sites, param);
		TimeGates(begin + 1, end, tau, sites, param);
	} else {
//		cout << "trotter 2 scheme" << endl;
//		double a1 = (3.-sqrt(3.))/6.0;		// more precise arrpoximation coefficients
//		double a2 = 1.0 - 2.0 * a1;
//		double b1 = 0.5;
//
//		double begin0 = begin; //this variable are needed to change operators A and B
//		double begin1 = begin + 1;
//		//Trotter gates from arxiv.org/abs/1901.04974
//		// Eq. (21)
//		cout << "Time evolutions " << endl;
//		TimeGates(begin0, end, a1 * tau, sites, param); //a1
//		TimeGates(begin1, end, b1 * tau, sites, param); //b1
//		TimeGates(begin0, end, a2 * tau, sites, param); //a2
//		TimeGates(begin1, end, b1 * tau, sites, param); //b1
//		TimeGates(begin0, end, a1 * tau, sites, param); //a1


		cout << "trotter 2 scheme" << endl;
		double begin0 = begin; //this variable are needed to change operators ABC
		double begin2 = begin + 1;
		double begin4 = begin + 2;
		//Trotter gates from arxiv.org/abs/1901.04974
		// Eq. (38),(47)

		cout << "Time evolutions " << endl;
		TimeGates(begin0, end, 0.5 * tau, sites, param); //A
		TimeGates(begin2, end, 0.5 * tau, sites, param); //B
		TimeGates(begin4, end, tau, sites, param); //C
		TimeGates(begin2, end, 0.5 * tau, sites, param); //B
		TimeGates(begin0, end, 0.5 * tau, sites, param); //A

	}
}

void TrotterExpXXZ::TimeGates(const int begin, const int end,
		const complex<double> tau, const SiteSet &sites,
		const ThreeSiteParam &param) {
	const int step = 3;
	const double J = param.val("J");
	const double Delta = param.val("Delta");
	const double J2 = param.val("J2");
	for (int j = begin; j < end - 1; j += step) {

		auto sp0sm1 = op(sites, "Sp", j) * op(sites, "Sm", j + 1) *
				op(sites, "Id", j + 2);
		auto sm0sp1 = op(sites, "Sm", j) * op(sites, "Sp", j + 1) *
				op(sites, "Id", j + 2);
		auto sz0sz1 = op(sites, "Sz", j) * op(sites, "Sz", j + 1) *
				op(sites, "Id", j + 2);
		auto sz0sz2 = op(sites, "Sz", j) * op(sites, "Id", j + 1) *
				op(sites, "Sz", j + 2);

		auto hh = J * 4 * 0.5 * sp0sm1;
		hh += J * 4 * 0.5 * sm0sp1;
		hh += Delta * 4 * sz0sz1;
		hh += J2 * 4 * sz0sz2;

		if(j == end - 2){
			hh += J * 4 * 0.5 * op(sites, "Id", j) *
					op(sites, "Sp", j + 1) * op(sites, "Sm", j + 2) ;
			hh += J * 4 * 0.5 * op(sites, "Id", j) *
					op(sites, "Sm", j + 1) * op(sites, "Sp", j + 2) ;
			hh += Delta * 4 * op(sites, "Id", j) *
					op(sites, "Sz", j + 1) * op(sites, "Sz", j + 2) ;
		}

		auto G = expHermitian(hh, tau);
		gates.emplace_back(j, move(G));
	}
}
void TrotterExpXXZ::Evolve(MPS &psi, const Args &args1) {
	const Args args("Truncate=", false);
	for (auto &gate : gates) {
		auto j = gate.i1;
		auto &G = gate.G;
		psi.position(j);
		auto WF = psi(j) * psi(j + 1) * psi(j + 2);
		WF = G * WF;
		WF /= norm(WF);
		WF.noPrime();
		{
			auto [Uj1, Vj1] = factor(WF,
					{ siteIndex(psi, j), leftLinkIndex(psi, j) }, args);
			auto indR = commonIndex(Uj1, Vj1);
			auto [Uj2, Vj2] = factor(Vj1, { siteIndex(psi, j + 1), indR }, args);
			psi.set(j, move(Uj1));
			psi.set(j + 1, move(Uj2));
			psi.set(j + 2, move(Vj2));

		}
	}
}





