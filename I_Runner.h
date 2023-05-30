#pragma once

#ifdef COMPILING_GA_DLL
#define DISCO2_GA_API __declspec(dllexport)
#else
#define DISCO2_GA_API __declspec(dllimport)
#endif

#include <vector>
#include "I_RunCase.h"

namespace DISCO2_GA {

/**
 Pure abstract class to define an access point for GA to utilize multiple simulators
 use in combination with I_RunCase and I_SigDataConvertor 
 to generate an eval value for each chromosome
*/
class DISCO2_GA_API I_Runner {
protected:

	std::vector<I_RunCase*> m_cases;

public:
	
	/**
	 Resets this I_Runner for a new NGA/SGA run
	*/
	virtual void reset() {};

	/**
	 Sets the cases that this will run during the multi-thread run
	*/
	virtual void setCaseEditors(std::vector<I_RunCase*>& caseList) {
		this->m_cases = caseList;
	}

	/**
	 Runs all the I_RunCase* cases
	*/
	virtual void runAll() = 0;

};

}	//namespace DISCO2_GA