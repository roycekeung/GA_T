#pragma once

#include <memory>

#include <Scenario.h>
#include <MultiThreadRunner.h>

#include "Job.h"

/**
A sample function on how to use MultiThreadRunner
for your main optimisation run function
JLo
*/
void sampleFunction() {
	// --- --- --- --- --- init convertor --- --- --- --- ---

	//add this to as a member variable to GA_Module, it should default construct during class construction anyways
	SigGpPhaseConvertor m_Convertor;

	//to init the thing use this function, pass the actual storage location of the respective things in param by reference
	//m_Convertor.init(std::shared_ptr<DISCO2_API::Scenario> scn, int sigSetUsing, 
	//		std::vector<int>& jctNodeIdSeq, std::vector<int>& phaseNumPerJct,
	//		std::vector< std::vector< std::tuple<int, int, int, bool, bool> > >& phaseInfo,
	//		std::vector< std::vector<std::pair< std::vector<int>, int>>>& crossPhaseMinGreen,
	//		std::vector< std::tuple<int, int, int, int, bool, bool> >& intersectionInfo);

	// --- --- --- --- --- using MultiThreadRunner --- --- --- --- --- 

	//from your param
	int sigSetUsing = 0;
	std::shared_ptr<DISCO2_API::Scenario> scn = nullptr;
	int numOfThreads = 1;

	//your job list
	std::vector<Job*> jobList{};

	//create the Runner, no need to create per generation
	DISCO2_API::MultiThreadRunner multiThreadRunner{ numOfThreads, scn }; 

	for (;;) {	//these can be looped in generations
		multiThreadRunner.setCaseEditors({ jobList.begin(), jobList.end() });
		multiThreadRunner.runAll();

		//to access all the results
		for (auto job : jobList) {
			job->getEvalValue();
		}
	}


}
