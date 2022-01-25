#pragma once

#include "itensor/all.h"
#include <vector>
#include <map>
#include <string>
#include <complex>

using namespace itensor;
using namespace std;


//------------------------------------------------------------------
//The function below translate numbers (etc.) into character strings
//the second parameter (optional) is the precision (digits)

//template<class T>
//string to_string(const T &t, unsigned int precision = 0) {
//	stringstream ss;
//	if (precision > 0) ss.precision(precision);
//	ss << t;
//	return ss.str();
//}

double char2double(char *a);

class Parameters: public map<string, double> { // class Parameters inherits all methods from container map<str,double>
public:
	double val(string var_name) const ;
	long longval(string var_name) const ;
	void PRint(ostream &o) const ;
	void ReadArguments(int argc, char *argv[]);
};

// class of PUBLIC parameters
class ThreeSiteParam: public Parameters {
public:
	ThreeSiteParam();
};


//I'm creating a Folded XXZ 3-site Hamiltonian
class ThreeSiteHamiltonian {
public:
	int dot;
	AutoMPO ampo;
	ThreeSiteHamiltonian(const SiteSet &sites, const ThreeSiteParam &param);
private:
	int N;
	void init(const ThreeSiteParam &param);
};


// A dual Folded XXZ 3-site Hamiltonian
class XXZ {
public:
	int dot;
	AutoMPO ampo;
	XXZ(const SiteSet &sites, const ThreeSiteParam &param);
private:
	int N;
	void init(const ThreeSiteParam &param);
};

//Folded XYZ 3-site Hamiltonian with integrability breaking term
class HamiltonianFoldedXYZ{
public:
	int dot;
	AutoMPO ampo;
	HamiltonianFoldedXYZ(const SiteSet &sites, const ThreeSiteParam &param);
private:
	int N;
	void init(const ThreeSiteParam &param);
};


//Trotter Gates for the time evolution of Folded_XXZ
class TrotterExp {
public:
	struct TGate {
		int i1 = 0;
		ITensor G;
		TGate(int i1_, ITensor G_)
		: i1(i1_)
		, G(G_) {
		}
	};
	TrotterExp(const SiteSet &sites, const ThreeSiteParam &param, const complex<double> tau);
	void initialize(const SiteSet &sites, const ThreeSiteParam &param, const complex<double> tau);
	void TimeGates(const int begin, const int end, const complex<double> tau,
			const SiteSet &sites, const ThreeSiteParam &param);
	void Evolve(MPS &psi, const Args &args) ;
private:
	vector<TGate> gates;
};


vector<MPO> XXZ_time_evol(const SiteSet &sites, const ThreeSiteParam &param);

//Trotter Gates for the time evolution of XXZ
class TrotterExpXXZ {
public:
	struct TGate {
		int i1 = 0;
		ITensor G;
		TGate(int i1_, ITensor G_)
		: i1(i1_)
		, G(G_) {
		}
	};
	TrotterExpXXZ(const SiteSet &sites, const ThreeSiteParam &param, const complex<double> tau);
	void initialize(const SiteSet &sites, const ThreeSiteParam &param, const complex<double> tau);
	void TimeGates(const int begin, const int end, const complex<double> tau,
			const SiteSet &sites, const ThreeSiteParam &param);
	void Evolve(MPS &psi, const Args &args) ;
private:
	vector<TGate> gates;
};


//Trotter Gates Folded XYZ with integrability breaking term
// https://arxiv.org/abs/2110.11322 in SIGMA repr Eq(1)
class TrotterExp_FoldedXYZ {
public:
	struct TGate {
		int i1 = 0;
		ITensor G;
		TGate(int i1_, ITensor G_)
		: i1(i1_)
		, G(G_) {
		}
	};
	TrotterExp_FoldedXYZ(const SiteSet &sites, const ThreeSiteParam &param, const complex<double> tau);
	void initialize(const SiteSet &sites, const ThreeSiteParam &param, const complex<double> tau);
	void TimeGates(const int begin, const int end, const complex<double> tau,
			const SiteSet &sites, const ThreeSiteParam &param);
	void Evolve(MPS &psi, const Args &args) ;
private:
	vector<TGate> gates;
};




