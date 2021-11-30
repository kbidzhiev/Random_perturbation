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
	const double Dhar = param.val("Dhar");
	const double PPX = sqrt(2.0)*param.val("PPX");
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

		if (PPX !=0){
			cout << "PPX term is included" << endl;
			auto I_0 = op(sites, "Id", j);
			auto Z_0 = op(sites, "Sz", j);
			auto I_1 = op(sites, "Id", j + 1);
			auto Z_1 = op(sites, "Sz", j + 1);
			auto Sx_2  = op(sites, "Sx", j + 2);

			hh += PPX * 0.5 * I_0 * I_1 * Sx_2;
			hh += -PPX * I_0 * Z_1 * Sx_2;
			hh += -PPX * Z_0 * I_1 * Sx_2;
			hh += PPX * 2.0 * Z_0 * Z_1 * Sx_2;
		}

		//Deepak Dhar term
		if (Dhar > 0) {
			cout << "Dhar term is included" << endl;
			hh += Dhar * 4 * 0.5 * op(sites, "Sz", j) * op(sites, "Id", j + 1)
					* op(sites, "Sz", j + 2);
			hh += -Dhar * 8 * 0.5 * op(sites, "Sz", j) * op(sites, "Sz", j + 1)
					* op(sites, "Sz", j + 2);
			hh += -Dhar * 1 * 0.5 * op(sites, "Id", j) * op(sites, "Id", j + 1)
					* op(sites, "Id", j + 2);
			hh += Dhar * 2 * 0.5 * op(sites, "Id", j) * op(sites, "Sz", j + 1)
					* op(sites, "Id", j + 2);
		}


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
			psi.set(j, Uj1);
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



//Trotter Gates for Exp_B
Exp_B::Exp_B(const SiteSet &sites, const ThreeSiteParam &param,
		const complex<double> tau) {
	initialize(sites, param, tau);
}
void Exp_B::initialize(const SiteSet &sites, const ThreeSiteParam &param,
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

void Exp_B::TimeGates(const int begin, const int end,
		const complex<double> tau, const SiteSet &sites,
		const ThreeSiteParam &param) {
	const int gate_range = 3;

	// 1/8 is a prefactor of exponent, 8 comes from SPin to Pauli
	// 1/2 from {SxSy - SySx} -> 0.5{SpSm-SmSp}
	const double coeff = 1.0; // 8 from spin matrices 8 from 4*2, so 8/8 == 1
	const double Delta_inverse = 1.0/param.val("Delta");


	for (int j = begin; j < end - 1; j += gate_range) {
		auto hh = Delta_inverse * coeff
				* op(sites, "Sx", j)
				* op(sites, "Sy", j + 1)
				* op(sites, "Sz", j + 2);

		hh += Delta_inverse * coeff
				* op(sites, "Sy", j)
				* op(sites, "Sx", j + 1)
				* op(sites, "Sz", j + 2);

		hh -= Delta_inverse * coeff
				* op(sites, "Sz", j)
				* op(sites, "Sx", j + 1)
				* op(sites, "Sy", j + 2);

		hh += Delta_inverse * coeff
				* op(sites, "Sz", j)
				* op(sites, "Sy", j + 1)
				* op(sites, "Sx", j + 2);

		auto G = expHermitian(hh, tau);
		gates.emplace_back(j, move(G));
	}
}
void Exp_B::Evolve(MPS &psi, const Args &args) {
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
			psi.set(j, Uj1);
			psi.set(j + 1, Uj2);
			psi.set(j + 2, Vj2);

		}
	}
}



vector<MPO> XP_time_evol(const SiteSet &sites, const ThreeSiteParam &param) {

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

	XP_Hamiltonian H_XP(sites, param);
	for (auto time_step : time_steps){
		MPO expH_i = toExpH(H_XP.ampo, time_step);
		Exp_H_vec.push_back(expH_i);
	}

	return Exp_H_vec;
}




//!!!!!!!!!!
//Trotter Gates PPK. ie I'm adding some perturbation to the dynamical Hamiltonian
TrotterExp_PPK::TrotterExp_PPK(const SiteSet &sites, const ThreeSiteParam &param,
		const complex<double> tau) {
	initialize(sites, param, tau);
}
void TrotterExp_PPK::initialize(const SiteSet &sites, const ThreeSiteParam &param,
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
		TimeGates(begin + 3, end, tau, sites, param);

	} else {
		cout << "trotter 2 scheme" << endl;
		cout << "its NOT implemented"<< endl;
		exit(1);

		/*
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
		 */
	}
}

void TrotterExp_PPK::TimeGates(const int begin, const int end,
		const complex<double> tau, const SiteSet &sites,
		const ThreeSiteParam &param) {
	const int gate_range = 4;
	const double J = param.val("J");
	const double PPK = sqrt(2.0)*param.val("PPK");
	//cout << "Gates starts from " << begin << endl;
	for (int j = begin; j < end - 2; j += gate_range) {
		//cout << "j = (" << j << "..., " << j + 3 <<  ")" << endl;

		auto I_0 = op(sites, "Id", j);
		auto Z_0 = op(sites, "Sz", j);
		auto Sp_0  = op(sites, "Sp", j);
		auto Sm_0  = op(sites, "Sm", j);
		auto I_1 = op(sites, "Id", j + 1);
		auto Z_1 = op(sites, "Sz", j + 1);
		auto Sp_2  = op(sites, "Sp", j + 2) * op(sites, "Id", j + 3);
		auto Sm_2  = op(sites, "Sm", j + 2) * op(sites, "Id", j + 3);
		//auto I_3   = op(sites, "Id", j + 3);

		auto hh = J * Sp_0 * I_1 * Sm_2 ;
		   hh +=  J * Sm_0 * I_1 * Sp_2 ;
		   hh += -J * 2 * Sp_0 * Z_1 * Sm_2 ;
		   hh += -J * 2 * Sm_0 * Z_1 * Sp_2 ;

		auto Kin_term1 = op(sites, "Sp", j + 2) * op(sites, "Sm", j + 3);
		auto Kin_term2 = op(sites, "Sm", j + 2) * op(sites, "Sp", j + 3);

		hh += PPK * 0.5 * I_0 * I_1 * Kin_term1; // (I)
		hh += PPK * 0.5 * I_0 * I_1 * Kin_term2; // (I)
		hh += -PPK * Z_0 * I_1 * Kin_term1; // (II)
		hh += -PPK * Z_0 * I_1 * Kin_term2; // (II)
		hh += -PPK * I_0 * Z_1 * Kin_term1; // (III)
		hh += -PPK * I_0 * Z_1 * Kin_term2; // (III)
		hh += PPK * 2.0 * Z_0 * Z_1 * Kin_term1; // (IV)
		hh += PPK * 2.0 * Z_0 * Z_1 * Kin_term2; // (IV)

		auto G = expHermitian(hh, tau);
		gates.emplace_back(j, move(G));
	}
}
void TrotterExp_PPK::Evolve(MPS &psi, const Args &args) {
	for (auto &gate : gates) {
		auto j = gate.i1;
		//cout << "j = " << j << endl;
		auto &G = gate.G;
		psi.position(j);
		auto WF = psi(j) * psi(j + 1) * psi(j + 2) * psi(j + 3);
		WF = G * WF;
		WF /= norm(WF);
		WF.noPrime();
		{
			auto [Uj1, Vj1] = factor(WF,
					{ siteIndex(psi, j), leftLinkIndex(psi, j) }, args);
			auto indR = commonIndex(Uj1, Vj1);
			auto [Uj2, Vj2] = factor(Vj1, { siteIndex(psi, j + 1), indR }, args);
			auto indR2 = commonIndex(Uj2, Vj2);
			auto [Uj3, Vj3] = factor(Vj2, { siteIndex(psi, j + 2), indR2 }, args);

			psi.set(j, Uj1);
			psi.set(j + 1, Uj2);
			psi.set(j + 2, Uj3);
			psi.set(j + 3, Vj3);

		}
		//cout << "GATES ARE APPLIED" << endl;
	}
}





