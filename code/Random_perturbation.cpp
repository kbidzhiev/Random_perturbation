#include "itensor/all.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
#include <sstream>
#include <complex>
#include <math.h>
#include <chrono>
#include <cmath>
#include <stdexcept>


#include <cstdlib>
#include "time_evolution.hpp"
#include "observables.hpp"
#include "quantum_measurement.hpp"

using namespace itensor;
using namespace std;
using namespace std::chrono;




//------------------------------------------------------------------
class MyDMRGObserver: public DMRGObserver { // vector<int> ~ DMRGObserver<ITensor>, where ITensor are tensor object type.
	double previous_energy;
	const double precision;
public:
	MyDMRGObserver(const MPS &psi, double prec = 1e-10) :
			DMRGObserver(psi), precision(prec) {
	}
	bool checkDone(const Args &args = Args::global()) { const
		double energy = args.getReal("Energy", 0);
		cout << "    Energy change:" << energy - previous_energy << endl;
		if (abs(energy - previous_energy) < precision) {
			cout << "   Energy has converged -> stop the DMRG.\n";
			return true;
		} else {
			previous_energy = energy;
			return false;
		}
	}
};

// main
int main(int argc, char *argv[]) {

	ThreeSiteParam param;
	param.ReadArguments(argc, argv); //Now param contains the parameters, default values or those provided on the command-line

	param.PRint(cout); // Print parameters
	cout.precision(15);

	const int N = param.longval("N");

	SpinHalf sites(N, { "ConserveQNs=", false }); // Creating a Hilbert space
	MPS psi, psi0;
	auto args = Args("Method=", "DensityMatrix", "Cutoff", param.val("trunc"),
			"MaxDim", param.longval("max_bond"), "Normalize", false); // for FitApplyMPO RENAME


// Preparing an initial state

	double energy_initial = 0;

	if (param.longval("GroundState") == 1) {
		// GS of the initial Hamiltonian
		cout << "initial state is GS" << endl;
		auto sweeps = Sweeps(999); //number of sweeps is 5
		sweeps.maxdim() = 10, 20, 50, 50, 100, 300, 4000;
		sweeps.cutoff() = 1E-10;

		psi = randomMPS(sites);

		cout << "Before DMRG" << endl;
		ThreeSiteHamiltonian Init_H(sites, param);
		auto H0 = toMPO(Init_H.ampo);


		cout << "Norm is = " << real(innerC(psi, psi)) << endl;
		energy_initial = real(innerC(psi, H0, psi)); //<psi|H0|psi>
		cout << "1. Initial energy=" << energy_initial << endl;

		MyDMRGObserver obs(psi, param.val("energy"));

		tie(energy_initial, psi) = dmrg(H0, psi, sweeps, obs, "Quiet");

		cout << "After DMRG" << endl;

		//int central_site = Init_H.dot;
		//SigmaXGate(central_site );
		//HadamarGate(central_site);
		psi.noPrime();

		//psi0 = psi; //we create two states. Psi for my time evolution, psi0 for standard one

	} else if (param.val("UDD") > 0) {
		cout << "initial state is  |UDD> " << endl;
		auto initState = InitState(sites);
		for (int i = 1; i <= N; ++i) {
			if (i % 3 == 0) {
				initState.set(i, "Up");
			} else {
				initState.set(i, "Dn");
			}
		}
		psi = MPS(initState);
		psi.noPrime();


		if(param.val("Perturb") != 0)
		{
			auto SigmaXGate = [&](int i) {
				auto ind = sites(i);
				auto indP = prime(sites(i));
				auto Op = ITensor(ind, indP);
				Op.set(ind(1), indP(1), 0);
				Op.set(ind(1), indP(2), 1);
				Op.set(ind(2), indP(1), 1);
				Op.set(ind(2), indP(2), 0);
				psi.setA(i, psi.A(i) * Op);
			};
			SigmaXGate(N / 2);
			psi.noPrime();
			cout << "spin is flipped" << endl;
		}
	} else if (param.val("UDDD") > 0) {
		cout << "initial state is  |UDDD> " << endl;
		auto initState = InitState(sites);
		for (int i = 1; i <= N; ++i) {
			if (i % 4 == 0) {
				initState.set(i, "Up");
			} else {
				initState.set(i, "Dn");
			}
		}
		psi = MPS(initState);
		psi.noPrime();

		if (param.val("Perturb") != 0) {
			auto SigmaXGate = [&](int i) {
				auto ind = sites(i);
				auto indP = prime(sites(i));
				auto Op = ITensor(ind, indP);
				Op.set(ind(1), indP(1), 0);
				Op.set(ind(1), indP(2), 1);
				Op.set(ind(2), indP(1), 1);
				Op.set(ind(2), indP(2), 0);
				psi.setA(i, psi.A(i) * Op);
			};
			SigmaXGate(N / 2);
			psi.noPrime();
			cout << "spin is flipped" << endl;
		}

	} else if ( param.val("RRR") > 0 ) {
		cout << "initial state is  |RRR> " << endl;
		auto initState = InitState(sites);
		for (int i = 1; i <= N; ++i){
				initState.set(i, "Up");
		}

		psi = MPS(initState);
		psi.noPrime();

		auto HadamarGate = [&](int i) {
			auto ind = sites(i);
			auto indP = prime(sites(i));
			auto Had = ITensor(ind, indP);
			Had.set(ind(1), indP(1), ISqrt2);
			Had.set(ind(1), indP(2), ISqrt2);
			Had.set(ind(2), indP(1), ISqrt2);
			Had.set(ind(2), indP(2), -ISqrt2);
			psi.setA(i, psi.A(i) * Had);
		};
		for (int i = 1; i <= N; ++i){
			HadamarGate(i);
		}
		psi.noPrime();

		if(param.val("Perturb") > 0) {
			RandomUnitary2sites(psi, sites, N/2, args);
		}
	} else if (param.val("Tilted") > 0) {
		cout << "initial state is  |RRR> " << endl;
		auto initState = InitState(sites);
		for (int i = 1; i <= N; ++i) initState.set(i, "Up");
		psi = MPS(initState);
		psi.noPrime();
		auto AlphaGate = [&](int i, const double alpha) {
			// U = exp[i alpha (n*s)]; |n|^2==1, s = {sx,sy,sz}
			// n = {0,1,0}
			auto ind = sites(i);
			auto indP = prime(sites(i));
			auto Op = ITensor(ind, indP);
			Op.set(ind(1), indP(1), cos(alpha));
			Op.set(ind(1), indP(2), sin(alpha));
			Op.set(ind(2), indP(1), -sin(alpha));
			Op.set(ind(2), indP(2), cos(alpha));
			psi.setA(i, psi.A(i) * Op);
		};
		for (int i = 1; i <= N; ++i) AlphaGate(i, param.val("Tilted"));

		auto SigmaXGate = [&](int i) {
			auto ind = sites(i);
			auto indP = prime(sites(i));
			auto Op = ITensor(ind, indP);
			Op.set(ind(1), indP(1), 0);
			Op.set(ind(1), indP(2), 1);
			Op.set(ind(2), indP(1), 1);
			Op.set(ind(2), indP(2), 0);
			psi.setA(i, psi.A(i) * Op);
		};
		if(param.val("Perturb") != 0) SigmaXGate(N / 2);
		psi.noPrime();

	} else if ( param.val("UUU") > 0 ) {
		cout << "initial state is  |uuuDuuu> " << endl;
		auto initState = InitState(sites);
		for (int i = 1; i <= N; ++i){
			if(i == N/2 + 1){
				initState.set(i, "Dn");
			} else {
				initState.set(i, "Up");
			}
		}
		psi = MPS(initState);
		psi.noPrime();


	} else if (param.val("Beta") > 0) {
		cout << "initial state is  |Exp[beta sum s_i*s_{i+1}]  > " << endl;
		auto initState = InitState(sites);
		for (int i = 1; i <= N; ++i)
			initState.set(i, "Up");
		psi = MPS(initState);
		psi.noPrime();


		double beta = param.val("Beta");
		double coshd = cosh(beta)/exp(beta);
		double sinhd = sinh(beta)/exp(beta);

		auto comInd = Index(2, "Link"); // define a common index

		for (int i = 1; i <= N; i++) { // run over the chain of size N
			auto s = sites(i);
			if (i == 1) {
				auto b = Index(2, "Link");
				auto A = ITensor(s, b);
				// here set the values of the matrix elements of A
				A.set(s(1), b(1), sqrt(coshd));
				A.set(s(1), b(2), sqrt(sinhd));
				A.set(s(2), b(1), sqrt(coshd));
				A.set(s(2), b(2), -sqrt(sinhd));
				psi.set(i, A);
				comInd = b; // set the common index to the right index of A
			}
			if (i > 1 && i < N) {
				auto a = comInd;
				auto b = Index(2, "Link");
				auto A = ITensor(s, a, b);

				A.set(s(1), a(1), b(1), coshd);
				A.set(s(1), a(1), b(2), sqrt(sinhd * coshd));
				A.set(s(1), a(2), b(1), sqrt(sinhd * coshd));
				A.set(s(1), a(2), b(2), sinhd);

				A.set(s(2), a(1), b(1), coshd);
				A.set(s(2), a(1), b(2), -sqrt(sinhd * coshd));
				A.set(s(2), a(2), b(1), -sqrt(sinhd * coshd));
				A.set(s(2), a(2), b(2), sinhd);

				psi.set(i, A);
				comInd = b; // set the common index to the right index of A
			}
			if (i == N) {
				auto a = comInd;
				auto A = ITensor(s, a);
				// here set the values of the matrix elements of A
				A.set(s(1), a(1), sqrt(coshd));
				A.set(s(1), a(2), sqrt(sinhd));
				A.set(s(2), a(1), sqrt(coshd));
				A.set(s(2), a(2), -sqrt(sinhd));
				psi.set(i, A);
			}


		} // end of for over sites




		if (param.val("Perturb") != 0){
			auto SigmaXGate = [&](int i) {
				auto ind = sites(i);
				auto indP = prime(sites(i));
				auto Op = ITensor(ind, indP);
				Op.set(ind(1), indP(1), 0);
				Op.set(ind(1), indP(2), 1);
				Op.set(ind(2), indP(1), 1);
				Op.set(ind(2), indP(2), 0);
				psi.setA(i, psi.A(i) * Op);
			};
			SigmaXGate(N / 2);
		}
		psi.noPrime();


	}else {
		cout << "Choose: GroundState, Neel, DomainWall,Impurity, Jammed = 1" << endl;
		return 1;
	}
	//cout << "PSI = " << psi << endl;
	psi0 = psi; //here I save the initial WF



//--------------------------------------------------------------

	//Hamiltonian for the dynamics
	//ThreeSiteHamiltonian Ham(sites, param);
	//auto H = toMPO(Ham.ampo);
	//if(param.val("Tilted") != 0){
	HamiltonianFoldedXYZ Ham(sites, param);
	auto H = toMPO(Ham.ampo);
	//}


	const int dot = Ham.dot;

	cout << "H constructed" << endl;

	energy_initial = real(innerC(psi, H, psi)); //<psi|H0|psi>
	cout << "2. Initial Energy psi =" << energy_initial << endl;

	double norm_initial = real(innerC(psi, psi)); //<psi|H0|psi>
	cout << "2. Initial NORM psi =" << norm_initial << endl;

// Output and observables
// _______
	ofstream ent, spec, entropy_profile, sz, sz_avrg, energy_beta,
				q1_profile, q2_profile, sx, sy, sxsysz,  sx_avrg, loschmidt,
				correlation1, correlation2, correlation3, sxsx, sztotal_strm; //here I'm defining output streams == files
	ios_base::openmode mode;
	mode = std::ofstream::out; //Erase previous file (if present)

	{
		sztotal_strm.open("Sz_total.dat", mode);
		sztotal_strm.precision(15);
		sztotal_strm << "#time \t  <Sz_total>  \t <Sz_total - sz_min> \t <Energy(t) - Energy0 >\n";
	}
	double dt = param.val("Loschmidt");
	if (dt != 0) { //Loschmidt echo
		loschmidt.open("Loschmidt.dat", mode);
		loschmidt.precision(15);
		loschmidt << "#time \t  re<psi(t)|psi(0)> \t re<psi(t)|psi(0)> \n";
	}
	//---------------------
	dt = param.val("Entropy");
	if (dt != 0) { //Entropy in the center of the chain
		ent.open("Entropy_center.dat", mode);
		ent.precision(15);
		ent << "#time \t Entropy(dot) \t BondDim(dot) \t MaxBondDim\n";
	}
	//---------------------
	dt = param.val("SVD_spec");
	if (dt > 0) { //SVD Spectrum on central bond
		spec.open("SVD_spec.dat", mode);
		spec.precision(15);
		spec << "#Position=" << dot << "\t<SVD_spectrum>\t\ttime\n";
	}

	//---------------------
	dt = param.val("EntropyProfile");
	if (dt > 0) { //Full entropy Profile
		entropy_profile.open("Entropy_profile.dat", mode);
		entropy_profile.precision(15);
		entropy_profile << "#Position=i-" << dot << setw(16) << "\t Entropy(i)"
				<< setw(16)
				<< "\t Entropy_sqrt \t Entropy_state1 \t time \t\t Bond.Dim(i)\n";
	}
	//---------------------
	dt = param.val("Sz");
	if (dt > 0) { //Full magnetization Profile
		sx.open("Sx_profile.dat", mode);
		sy.open("Sy_profile.dat", mode);
		sz.open("Sz_profile.dat", mode);
		sxsysz.open("SxSySz_profile.dat", mode);
		sxsx.open("SxSx_r_correlation.dat", mode);

		correlation1.open("Correlation1.dat", mode);
		correlation2.open("Correlation2.dat", mode);
		correlation3.open("Correlation3.dat", mode);

		sz_avrg.open("Sz_average_profile.dat", mode);
		sx_avrg.open("Sx_average_profile.dat", mode);

		sx.precision(15);
		sy.precision(15);
		sz.precision(15);
		sxsysz.precision(15);
		sxsx.precision(15);


		sz_avrg.precision(15);
		sx_avrg.precision(15);

		sx << "#Position=i-"
				<< "\t <Sx_i>\t"
				<< "\t <Sx_i Sx_{i+1}>\t"
				<< "\t <Sx_i Sx_{i+2}>\t"
				<< "\t <Sx_i Sx_{i+3}>\t"
				<< "\t <Sx_i Sx_{i+4}>\t"
				<< "\t <Sx_i Sz_{i+1}>\t"
				<< "\t\ttime\n";
		sy << "#Position=i-"
				<< "\t <Sy_i>\t"
				<< "\t <Sy_i Sy_{i+1}>\t"
				<< "\t <Sy_i Sy_{i+2}>\t"
				<< "\t <Sy_i Sy_{i+3}>\t"
				<< "\t <Sy_i Sy_{i+4}>\t"
				<< "\t <Sy_i Sz_{i+1}>\t"
				<< "\t\ttime\n";
		sz << "#Position=i-"
				<< "\t <Sz_i>\t"
				<< "\t <Sz_i Sz_{i+1}>\t"
				<< "\t <Sz_i Sz_{i+2}>\t"
				<< "\t <Sz_i Sz_{i+3}>\t"
				<< "\t <Sz_i Sz_{i+4}>\t"
				<< "\t (-1)^i<Sz_i>\t"
				<< "\t\ttime\n";

		sxsysz << "#Position=i-"
				<< "\t <Sx_i>\t"
				<< "\t <Sy_i>\t"
				<< "\t <Sz_i>\t"
				<< "\t\ttime\n";

		sxsx	<< "#Position=i-"
				<< "\t <SxSx>\t"
				<< "\t RE<SxSy>\t"
				<< "\t IM<SxSy>\t"
				<< "\t\ttime\n";

		correlation1<< "#Position=i-"
				<< "\t <Sx_i Sx_{i+1}>\t"	//2
				<< "\t <Sx_i Sy_{i+1}>\t"	//3
				<< "\t <Sx_i Sz_{i+1}>\t"	//4
				<< "\t <Sy_i Sx_{i+1}>\t"	//5
				<< "\t <Sy_i Sy_{i+1}>\t"	//6
				<< "\t <Sy_i Sz_{i+1}>\t"	//7
				<< "\t <Sz_i Sx_{i+1}>\t"	//8
				<< "\t <Sz_i Sy_{i+1}>\t"	//9
				<< "\t <Sz_i Sz_{i+1}>\t"	//10
				<< "\t\ttime\n";

		correlation2<< "#Position=i-"
				<< "\t <Sx_i Sx_{i+2}>\t" 	//2
				<< "\t <Sx_i Sy_{i+2}>\t"	//3
				<< "\t <Sx_i Sz_{i+2}>\t"	//4
				<< "\t <Sy_i Sx_{i+2}>\t"	//5
				<< "\t <Sy_i Sy_{i+2}>\t"	//6
				<< "\t <Sy_i Sz_{i+2}>\t"	//7
				<< "\t <Sz_i Sx_{i+2}>\t"	//8
				<< "\t <Sz_i Sy_{i+2}>\t"	//9
				<< "\t <Sz_i Sz_{i+2}>\t"	//10
				<< "\t\ttime\n";

		correlation3<< "#Position=i-"
				<< "\t <Sx_i Sx_{i+3}>\t" 	//2
				<< "\t <Sx_i Sy_{i+3}>\t"	//3
				<< "\t <Sx_i Sz_{i+3}>\t"	//4
				<< "\t <Sy_i Sx_{i+3}>\t"	//5
				<< "\t <Sy_i Sy_{i+3}>\t"	//6
				<< "\t <Sy_i Sz_{i+3}>\t"	//7
				<< "\t <Sz_i Sx_{i+3}>\t"	//8
				<< "\t <Sz_i Sy_{i+3}>\t"	//9
				<< "\t <Sz_i Sz_{i+3}>\t"	//10
				<< "\t\ttime\n";

		sz_avrg << "#Position=i-"
				<< "\t 0.5*<Sz_i + next site>\t"
				<< "\t 0.5*<Sz_i Sz_{i+1} + next site>\t"
				<< "\t 0.5*<Sz_i Sz_{i+2} + next site>\t"
				<< "\t 0.5*<Sz_i Sz_{i+3} + next site>\t"
				<< "\t 0.5*<Sz_i Sz_{i+4} + next site>\t"
				<< "\t\ttime\n";

		sx_avrg << "#Position=i-"
				<< "\t 0.5*<Sx_i + next site>\t"
				<< "\t 0.5*<Sx_i Sx_{i+1} + next site>\t"
				<< "\t 0.5*<Sx_i Sx_{i+2} + next site>\t"
				<< "\t 0.5*<Sx_i Sx_{i+3} + next site>\t"
				<< "\t 0.5*<Sx_i Sx_{i+4} + next site>\t"
				<< "\t\ttime\n";
	}
	//---------------------
	dt = param.val("Q1Profile");
	if (dt > 0) { //Energy Profile and Current profile
		q1_profile.open("Q1_profile.dat", mode);
		q1_profile.precision(15);
		q1_profile << "#Position=i-" << "\t<Q1minus_i>\t" << dot
				<< "\t\ttime(or beta)\n";

	}
	//---------------------
	dt = param.val("Q2Profile");
	if (dt > 0) { //Full Q2 Profile
		q2_profile.open("Q2_profile.dat", mode);
		q2_profile.precision(15);
		q2_profile << "#Position=i-" << dot << setw(16) << "\t Entropy(i)"
				<< setw(16) << "\t Q2plus \t Q2minus \t time \t \n";
	}

	//exp(Ham) for the dynamics


	double tau = param.val("tau");

	const long int n_steps = param.val("T") / tau;
	TrotterExp expH_Folded_XXZ(sites, param, -Cplx_i * tau);
	TrotterExpXXZ expH_XXZ(sites, param, -Cplx_i * tau);
	TrotterExp_FoldedXYZ expH_XYZ(sites, param, -Cplx_i * tau);
	vector<MPO> XXZ_time_evol_vec = XXZ_time_evol(sites, param);



// Time evolution
	double sz_total_initial = 9999;
	double sz_tot = 999; // will be used in part with magnetization and then output to cout;
	double sx_tot = 999; // will be used in part with magnetization and then output to cout;
	for (int n = 0; n <= n_steps ; ++n) {

		const double time = n * tau; //+param.val("time_shift");
		cout << "Time step #" << n << "/" << n_steps << "\ttime=" << time
				<< endl;

		vector<double> Myspec; //vector which will be the SVD spectrum

		//  ----- Compute observables ---
		//  ----- Entropy between sites i and i+1
		if (param.val("Entropy") != 0) {
			double entr = Entropy(psi, dot, Myspec, 1);
			ent << time << "\t" << setw(16) << setfill('0') << entr << "\t"
					<< BondDim(psi, dot) << "\t" << maxLinkDim(psi) << endl;

			if (param.val("SVD_spec") > 0) {
				if (n % int(param.val("SVD_spec") / param.val("tau")) == 0) {
					spec << "\"t=" << time << "\"" << endl;
					int si = Myspec.size();
					for (int i = 0; i < si; i++) {
						spec << i + 1 << "\t" << Myspec[i] << "\t"
								<< sqrt(Myspec[i]) << "\t" << time << endl;
					}
					if (n < n_steps)
						spec << endl << endl;
				}
			}
		}


		// ------- Loschmidt echo -----
		if (param.val("Loschmidt") > 0)
			if (n % int(param.val("Loschmidt") / param.val("tau")) == 0) {
				complex<double> echo = innerC(psi0,psi);
				loschmidt << time << "\t" << setw(16) << setfill('0')
						<< real(echo) << "\t" << imag(echo)	<< endl;
			}

		// ------- entanglement Entropy Profile -----
		if (param.val("EntropyProfile") > 0)
			if (n % int(param.val("EntropyProfile") / param.val("tau")) == 0) {
				entropy_profile << "\"t=" << time << "\"" << endl;
				for (int i = 1; i < N; i++) {
					double entr_std = Entropy(psi, i, Myspec, 1); // p log p
					entropy_profile << i + 0.5 - dot << "\t" << setw(16) << setfill('0')
							<< entr_std << "\t" << setw(16) << setfill('0')
							<< BondDim(psi, i) << "\t" << time << endl;
				}
				if (n < n_steps)
					entropy_profile << endl << endl;
			}
		// ------- Sz profile -------

		if (param.val("Sz") > 0) {
			if (n % int(param.val("Sz") / param.val("tau")) == 0) {
				sx << "\"t=" << time << "\"" << endl;
				sy << "\"t=" << time << "\"" << endl;
				sz << "\"t=" << time << "\"" << endl;
				sxsx << "\"t=" << time << "\"" << endl;

				sxsysz << "\"t=" << time << "\"" << endl;

				correlation1 << "\"t=" << time << "\"" << endl;
				correlation2 << "\"t=" << time << "\"" << endl;
				correlation3 << "\"t=" << time << "\"" << endl;

				sx_avrg << "\"t=" << time << "\"" << endl;
				sz_avrg << "\"t=" << time << "\"" << endl;

				double sz_left = 0, sz_right = 0, sz_dot = 0;
				sz_tot = 0;
				sx_tot = 0;

				double sx_odd = 0,  sz_odd = 0;
				double sxsx1_odd = 0, sxsx2_odd = 0, sxsx3_odd = 0, sxsx4_odd = 0;
				double szsz1_odd = 0, szsz2_odd = 0, szsz3_odd = 0, szsz4_odd = 0;

				double sxsz_odd  = 0; //sysz_odd  = 0;

				for (int i = 1; i <= N; i ++) {
					double s = Sz(psi, sites, i);
					double sx1 = Sx(psi, sites, i);
					double sy1 = Sy(psi, sites, i);
					double sxsx4 = 0;
					double sysy4 = 0;
					double szsz4 = 0;

					// some correlations have form C(i,i+2)
					// if i > N-2 I define correlation to be 0
					double sxsx1 = 0, sxsy1 = 0, sxsz1 = 0;
					double sysx1 = 0, sysy1 = 0, sysz1 = 0;
					double szsx1 = 0, szsy1 = 0, szsz1 = 0;
					double sxsx2 = 0, sxsy2 = 0, sxsz2 = 0;
					double sysx2 = 0, sysy2 = 0, sysz2 = 0;
					double szsx2 = 0, szsy2 = 0, szsz2 = 0;

					double sxsx3 = 0, sxsy3 = 0, sxsz3 = 0;
					double sysx3 = 0, sysy3 = 0, sysz3 = 0;
					double szsx3 = 0, szsy3 = 0, szsz3 = 0;

					complex<double> kss = 0;

					if (i <= N - 1) {
						sxsx1 = real(Correlation(psi, sites, "Sx", "Sx", i, i+1));
						sxsy1 = real(Correlation(psi, sites, "Sx", "Sy", i, i+1));
						sxsz1 = real(Correlation(psi, sites, "Sx", "Sz", i, i+1));

						sysx1 = real(Correlation(psi, sites, "Sy", "Sx", i, i+1));
						sysy1 = real(Correlation(psi, sites, "Sy", "Sy", i, i+1));
						sysz1 = real(Correlation(psi, sites, "Sy", "Sz", i, i+1));

						szsx1 = real(Correlation(psi, sites, "Sz", "Sx", i, i+1));
						szsy1 = real(Correlation(psi, sites, "Sz", "Sy", i, i+1));
						szsz1 = real(Correlation(psi, sites, "Sz", "Sz", i, i+1));
					}
					if (i <= N - 2) {
						sxsx2 = real(Correlation(psi, sites, "Sx", "Sx", i, i+2));
						sxsy2 = real(Correlation(psi, sites, "Sx", "Sy", i, i+2));
						sxsz2 = real(Correlation(psi, sites, "Sx", "Sz", i, i+2));

						sysx2 = real(Correlation(psi, sites, "Sy", "Sx", i, i+2));
						sysy2 = real(Correlation(psi, sites, "Sy", "Sy", i, i+2));
						sysz2 = real(Correlation(psi, sites, "Sy", "Sz", i, i+2));

						szsx2 = real(Correlation(psi, sites, "Sz", "Sx", i, i+2));
						szsy2 = real(Correlation(psi, sites, "Sz", "Sy", i, i+2));
						szsz2 = real(Correlation(psi, sites, "Sz", "Sz", i, i+2));

						kss = KSS(psi, sites, i);
					}
					const int distance = 3;
					if (i <= N - distance) {

						sxsx3 = real(Correlation(psi, sites, "Sx", "Sx", i, i+distance));
						sxsy3 = real(Correlation(psi, sites, "Sx", "Sy", i, i+distance));
						sxsz3 = real(Correlation(psi, sites, "Sx", "Sz", i, i+distance));
						sysx3 = real(Correlation(psi, sites, "Sy", "Sx", i, i+distance));
						sysy3 = real(Correlation(psi, sites, "Sy", "Sy", i, i+distance));
						sysz3 = real(Correlation(psi, sites, "Sy", "Sz", i, i+distance));
						szsx3 = real(Correlation(psi, sites, "Sz", "Sx", i, i+distance));
						szsy3 = real(Correlation(psi, sites, "Sz", "Sy", i, i+distance));
						szsz3 = real(Correlation(psi, sites, "Sz", "Sz", i, i+distance));

					}
					if (i <= N - 4) {
						sxsx4 = real(Correlation(psi, sites, "Sx", "Sx", i, i+4));
						sysy4 = real(Correlation(psi, sites, "Sy", "Sy", i, i+4));
						szsz4 = real(Correlation(psi, sites, "Sz", "Sz", i, i+4));
					}
					if (i > N / 2) {
						complex<double> sxsx_r = 0;
						complex<double> sxsy_r = 0;
						double szsz_r =0;

						sxsx_r = 	Correlation(psi, sites, "Sx", "Sx", N/2, i);
						sxsy_r = 	Correlation(psi, sites, "Sx", "Sy", N/2, i);
						szsz_r = real(Correlation(psi, sites, "Sz", "Sz", N/2, i));

						sxsx << i - dot << "\t"
								<< sxsx_r.real() << "\t"
								<< sxsx_r.imag() << "\t"
								<< sxsy_r.real() << "\t"
								<< sxsy_r.imag() << "\t"
								<< szsz_r << "\t"
								<< time	<< endl;
					}


					sz_tot += s;
					sx_tot += sx1;
					if (i < dot)
						sz_left += s;
					if (i > dot)
						sz_right += s;
					if (i == dot)
						sz_dot += s;
					sz << i - dot  << "\t"			//column ID below
						<< s << "\t"				//2
						<< szsz1 << "\t"			//3
						<< szsz2 << "\t"			//4
						<< szsz3 << "\t"			//5
						<< szsz4 << "\t"			//6
						<< pow(-1, i ) * s << "\t"	//7
						 << "\t"<< time << endl;
					sx << i - dot  << "\t"			//column ID below
						<< sx1 << "\t"				//2
						<< sxsx1 << "\t"			//3
						<< sxsx2 << "\t"			//4
						<< sxsx3 << "\t"			//5
						<< sxsx4 << "\t"			//6
						<< sxsz1 << "\t"				//7
						<< time << endl;
					sy << i - dot  << "\t"			//column ID below
						<< sy1 << "\t"				//2
						<< sysy1 << "\t"			//3
						<< sysy2 << "\t"			//4
						<< sysy3 << "\t"			//5
						<< sysy4 << "\t"			//6
						<< sysz1 << "\t"			//7
						<< time << endl;
					sxsysz << i - dot  << "\t"		//column ID below
						<< sx1 << "\t"				//2
						<< sy1 << "\t"				//3
						<< s << "\t"				//4
						<< real(kss) << "\t"		//5
						<< time << endl;

					correlation1 << i - dot  << "\t"	//column ID below
						<< sxsx1 << "\t"				//2
						<< sxsy1 << "\t"				//3
						<< sxsz1 << "\t"				//4
						<< sysx1 << "\t"				//5
						<< sysy1 << "\t"				//6
						<< sysz1 << "\t"				//7
						<< szsx1 << "\t"				//8
						<< szsy1 << "\t"				//9
						<< szsz1 << "\t"				//10
						<< time << endl;

					correlation2 << i - dot  << "\t"	//column ID below
						<< sxsx2 << "\t"				//2
						<< sxsy2 << "\t"				//3
						<< sxsz2 << "\t"				//4
						<< sysx2 << "\t"				//5
						<< sysy2 << "\t"				//6
						<< sysz2 << "\t"				//7
						<< szsx2 << "\t"				//8
						<< szsy2 << "\t"				//9
						<< szsz2 << "\t"				//10
						<< time << endl;

					correlation3 << i - dot  << "\t"	//column ID below
						<< sxsx3 << "\t"				//2
						<< sxsy3 << "\t"				//3
						<< sxsz3 << "\t"				//4
						<< sysx3 << "\t"				//5
						<< sysy3 << "\t"				//6
						<< sysz3 << "\t"				//7
						<< szsx3 << "\t"				//8
						<< szsy3 << "\t"				//9
						<< szsz3 << "\t"				//10
						<< time << endl;


					if ( i % 2 == 1) { //odd site
						sz_odd = s;
						sx_odd = sx1;
						szsz1_odd = szsz1;
						szsz2_odd = szsz2;
						szsz3_odd = szsz3;
						szsz4_odd = szsz4;

						sxsx1_odd = sxsx1;
						sxsx2_odd = sxsx2;
						sxsx3_odd = sxsx3;
						sxsx4_odd = sxsx4;
						sxsz_odd = sxsz1;
					} else {
						sz_avrg << i - dot + 1 << "\t"
								<< 0.5 * (s + sz_odd) << "\t"
								<< 0.5 * (szsz1 + szsz1_odd) << "\t"
								<< 0.5 * (szsz2 + szsz2_odd) << "\t"
								<< 0.5 * (szsz3 + szsz3_odd) << "\t"
								<< 0.5 * (szsz4 + szsz4_odd) << "\t"
								<< time << endl;
						sx_avrg << i - dot + 1 << "\t"
								<< 0.5 * (sx1 + sx_odd) << "\t"
								<< 0.5 * (sxsx1 + sxsx1_odd) << "\t"
								<< 0.5 * (sxsx2 + sxsx2_odd) << "\t"
								<< 0.5 * (sxsx3 + sxsx3_odd) << "\t"
								<< 0.5 * (sxsx4 + sxsx4_odd) << "\t"
								<< 0.5 * (sxsz1  + sxsz_odd) << "\t"
								<< time << endl;
					}
				}

				{ //I need this part to separate time steps in *.dat files (for gnuplot)
					sx << "\n\n";
					sx_avrg << "\n\n";
					sy << "\n\n";
					sz << "\n\n";
					sz_avrg << "\n\n";
					sxsysz << "\n\n";
					correlation1 << "\n\n";
					correlation2 << "\n\n";
					correlation3 << "\n\n";
					sxsx << "\n\n";
				}
				//sz_tot += Sz(psi, sites, N-1) + Sz(psi, sites, N);
				if (n == 0) {
					sz_total_initial = sz_tot;
				}
				cout << "\n<Sz_left>=" << sz_left << "\t" << "<Sz_right>="
						<< sz_right << "\t" << "<Sz_DOT>=" << sz_dot << "\t"
						<< "<Sz_tot>=" << sz_tot
						<<endl;
			}
		}
		// ------- Energy Profile -------
		if (param.val("Q1Profile") > 0 ) {
			if (n % int(param.val("Q1Profile") / param.val("tau")) == 0) {
				q1_profile << "\"t=" << time << "\"" << endl;
				for (int i = 1; i <= N - 3; i++) {
					const complex<double> q1 = Q1(psi, sites, i);
					const double en = real(q1);
					const double q1minus = imag(q1);
					q1_profile << i - dot
							<< "\t" << en
							<< "\t" << q1minus
							<< "\t" << time << endl;
				}
				q1_profile << "\n\n"; //I need this part to separate time steps in *.dat files (for gnuplot)
			}
		}
		// ------- Q2 Profile -------
		if (param.val("Q2Profile") > 0 ) {
			if (n % int(param.val("Q2Profile") / param.val("tau")) == 0) {
				q2_profile << "\"t=" << time << "\"" << endl;
				for (int i = 1; i <= N - 5; i++) {
					const complex<double> q2 = Q2(psi, sites, i);
					const double q2_plus =  real(q2);
					const double q2_minus = imag(q2);
					q2_profile << i - dot << "\t" << q2_plus << "\t"
							<< q2_minus << "\t" << time << endl;
				}
				q2_profile << "\n\n"; //I need this part to separate time steps in *.dat files (for gnuplot)
			}
		}

		if (n < n_steps) {
			//MPS psi_temp = psi;
			cout << "Time evol" << endl;
			if(param.val("XXZ") > 0 ){
				if(param.val("TrotterOrderXXZ") > 0){
					for(auto & expH_XXZ : XXZ_time_evol_vec){
						psi = applyMPO(expH_XXZ, psi, args);
					}
					psi.noPrime();
					cout << "W^{II} operator XXZ" << endl;
				} else{
					expH_XXZ.Evolve(psi, args);
					cout << "TEBD XXZ" << endl;
				}
			} else if (param.val("Tilted") > 0
						|| param.val("UDD") > 0
						|| param.val("UDDD") > 0
						|| param.val("UUU") > 0
						|| param.val("Beta") > 0) {
				expH_XYZ.Evolve(psi, args);
				cout << "Folded XYZ from Maurizios' paper" << endl;
			} else {
				expH_Folded_XXZ.Evolve(psi, args);
				cout << "Folded XXZ" << endl;
			}
			psi.orthogonalize(args);

			double sz_TOT = 0;
			for (int i = 1; i <= N; i++) {
				sz_TOT += Sz(psi, sites, i);
			}
			double sz_dif = sz_TOT - sz_total_initial;

			double energy = real(innerC(psi, H, psi));
			double energy_diff = energy - energy_initial;
			cout << "max bond dim = " << maxLinkDim(psi) << endl;
			//cout << "Norm = " << real(innerC(psi, psi)) << endl;
			cout << "Energy = " << energy << endl;
			cout << "E(t) - E(0) = " << energy_diff << endl;
			cout << "Sz(t) - Sz(0) = " << sz_dif << endl;

			sztotal_strm << time << "\t" << setw(16) //<< setfill('0')
					<< sz_TOT << "\t" << sz_dif << "\t" << energy_diff << endl;
		}
	}
	cout << "\nTime evolution complete.\n";
	cout << "Done !" << endl;
	return 0;
}

