#pragma once

#ifdef COMPILING_GA_DLL
#define DISCO2_GA_API __declspec(dllexport)
#else
#define DISCO2_GA_API __declspec(dllimport)
#endif

#include <unordered_map>
#include "GA_Structs.h"

namespace DISCO2_GA {

class Chromosome;

/**
 An abstract class to define a customisable sig data convertor
 to access/initiate the internal GA reference data
 and to interpret the Chromosome
 a pointer to this class (same instance) will be provided to all I_RunCase
*/
class DISCO2_GA_API I_SigDataConvertor {
protected:
	Chromosome* t_bestChromosome = nullptr;

public:

	/**
	 Function for inputting reference data into GA
	 -> see struct defs in GA_Structs.h
	*/
	virtual void init(std::unordered_map<int, intersec_inf>& intersecs_inf) = 0;

	/**
	 GA internal function for GA module to automatically insert the best chromo into convertor
	  intended for easier output/translation to elsewhere
	*/
	void setTChromosome(Chromosome* tChromo) {
		this->t_bestChromosome = tChromo;
	}

};

}	//namespace DISCO2_GA