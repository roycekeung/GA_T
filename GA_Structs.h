#pragma once

#ifdef COMPILING_GA_DLL
#define DISCO2_GA_API __declspec(dllexport)
#else
#define DISCO2_GA_API __declspec(dllimport)
#endif

#include <vector>
#include <utility>

namespace DISCO2_GA {

/**
 1 phase information in 1 intersection including:
 changeable, min_green, intergreen , cross phase (eg."0100110" or int 45)
*/
struct DISCO2_GA_API phase_inf {
	// a start up value of a phase duration in an intersection
	int initial_phase_duration = 0;

	// a determination whether this phase duration could be cange or not
	bool phase_changeable_by_value = true;   // if false, fixed by exact value
	bool phase_changeable_by_ratio = true;   // if false, fixed by ratio

	// minimum phase duration eg. 15s for each phase duration in 1 intersection
	int min_green = 0;

	// intergreen for each phase duration in 1 intersection
	int intergreen = 0;

	// the range of values phase duration allowed to be changed from the initial guess
	// limits the lower bound of the phase green time
	// minGreen will always be respected
	// input must be positive number, -1 -> unlimited change
	int changableInterval = -1;
};


/**
 1 intersection information including:
 typeOfplan, offset, cycle, changeable,
 min_cycle, max_cycle, ptr collections of all phase information, all in binary and decimal (eg."0100110" or int 45)
*/
struct DISCO2_GA_API intersec_inf {
	// type of timing plan by each intersection
	enum class GAtype_SigCtrl { FGFC, VGFC, VGVC, other };
	GAtype_SigCtrl typeOfplan_per_intersec;

	// how many number of cycle to fullfil the total simulation time under the VGVC or VGFC  typeOfplan
	int no_Cycle_intotal = 1;

	// a start up value of an offset in an intersection
	int initial_Offset;

	// a determination whether this offset could be cHangeD or not
	bool Offset_changeable;

	// a start up value of a cycle in an intersection
	int initial_Cycle;

	// a determination whether this cycle could be change or not
	bool Cycle_changeable;

	// maximum Cycle time eg. 120s usually in hk
	int max_Cycle;

	// minimum Cycle time per intersection, eg.  all min green + all intergreen in a intersection
	//(but now's mannually input, so it could be greater than actual  all min green + all intergreen )
	int min_Cycle;

	// the range of values phase duration allowed to be changed from the initial guess
	int changableInterval = -1;

	// collectinos of phase_duration information 
	std::vector<phase_inf> phase_duration_inf_collections;

	// compute this after the inf of phase table have been intactly described;
	// ideally ,  min_Cycle  = all min green(sum_PhaseMinG) + all intergreen (sum_PhaseInterG),   
	//	but they could be diff due to these two variable is cal later
	// they r used in the gen_random_timeplan() and the trans func 
	int sum_PhaseMinG;
	int sum_PhaseInterG;

	//for the cross phase setting (if the phases in side this intersection dont have crossphase, then tempPhasesMinGreen would be empty )
	//  std::pair first : list of phases index, second : crossphase min green 
	//	=====> this is for the crossphase only, the vector size depends on how many crossphase u have
	std::vector<std::pair<std::vector<int>, int>> PhaseidList_CrossPhasesMinGreen;


	// constructor of intersec_inf
	intersec_inf(GAtype_SigCtrl typeOfplan_per_intersec,
		int initial_Offset,
		bool Offset_changeable,
		int initial_Cycle,
		bool Cycle_changeable,
		int max_Cycle,
		int min_Cycle,
		std::vector<phase_inf> phase_duration_inf_collections,
		std::vector<std::pair<std::vector<int>, int>> PhaseidList_CrossPhasesMinGreen,
		int changableInterval = -1) :

		typeOfplan_per_intersec(typeOfplan_per_intersec),
		initial_Offset(initial_Offset),
		Offset_changeable(Offset_changeable),
		initial_Cycle(initial_Cycle),
		Cycle_changeable(Cycle_changeable),
		max_Cycle(max_Cycle),
		min_Cycle(min_Cycle),
		phase_duration_inf_collections(phase_duration_inf_collections),
		PhaseidList_CrossPhasesMinGreen(PhaseidList_CrossPhasesMinGreen),
		changableInterval(changableInterval) {

		sum_PhaseMinG = 0;
		sum_PhaseInterG = 0;

		for (phase_inf& phase : phase_duration_inf_collections) {
			sum_PhaseMinG += phase.min_green;
			sum_PhaseInterG += phase.intergreen;
		}
	}
};



}