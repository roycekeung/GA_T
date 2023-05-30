#pragma once

#ifdef COMPILING_GA_DLL
#define DISCO2_GA_API __declspec(dllexport)
#else
#define DISCO2_GA_API __declspec(dllimport)
#endif

#include "Chromosome.h"
#include "I_SigDataConvertor.h"

namespace DISCO2_GA {

/**
 An abstract class to define an access point for GA to link its Chromosome with I_SigDataConvertor
 such that each chromosome will generate an eval value for assessing the quality of chromosome
 use in combination with I_Runner to provide customisation of simulators
*/
class DISCO2_GA_API I_RunCase {
protected:

	double m_evalValue = -1;

	Chromosome m_chromo;

	const I_SigDataConvertor* ref_convertor;

public:

	// --- --- --- --- --- Constructor Destructor --- --- --- --- ---

	I_RunCase(Chromosome&& chromo, const I_SigDataConvertor* convertor) :
		m_chromo(chromo), ref_convertor(convertor) {}

	virtual ~I_RunCase() {}

	// --- --- --- --- --- Getters --- --- --- --- ---

	// get the Chromosome_obj out
	Chromosome& get_Chromosome_obj() { return this->m_chromo; }

	/**
	 Used for evaluating the chromosome
	*/
	double getEvalValue() { return this->m_evalValue; }

};

/**
 Factory of I_RunCase such that I_RunCase can be created from within GA
*/
class DISCO2_GA_API I_RunCase_Factory {
public:

	virtual I_RunCase* genRunCase(Chromosome&& chromo, const I_SigDataConvertor* convertor) = 0;

};

}	//namespace DISCO2_GA