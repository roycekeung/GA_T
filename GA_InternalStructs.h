#pragma once

#include <set>
#include <vector>
#include <unordered_map>
#include <limits>

namespace DISCO2_GA {

class I_RunCase;

/**
 one the intersection involve in the intersection group setting,
 then it would have it own intersecs_gps_inf struct
*/
struct intersecs_gps_inf {
	// the base coordinated intersection id, the cycle time is synchronize with this
	int coordinated_intersec_id = -1;

	// relationship list of all intersection id that r in the group
	std::vector <int> relation_group_id = {};

	// used for the worst case scenario, in the decoding translation part in oder to get the same cycle time value
	// min cycle_max among all relation_group_id intersections
	int minimum_maxCycle = std::numeric_limits<int>::max();
	// max cycle_min among all relation_group_id intersections
	int maximum_minCycle = 0;

};

struct input_param {
	// --- --- --- --- --- type of GA strategy --- --- --- --- --- //
	enum class GAtype { NGA, SGA };
	GAtype GA_type = GAtype::NGA;   // by default

	// --- --- --- --- --- compulsory input --- --- --- --- --- //
	int debugPrintOutLevel = 1;

	int trial_count_limit = 3;
	double Crossover_rate = 0.5;		// percentage of the bit-wise switch of parents take plpace
	double crossover_multiplier = 1;	// no_element * crossover_multiplier = the number of crossover point
	double Mutation_rate = 0.005;		// percentage of the bit-wise switch between 0,1 take plpace
	double Power_factor = 5;			// amplification for the fitness function
	int Populcation_size = 100;			// total number of chromosome in each generation run
	double Elite_rate = 0.08;			// percentage to total number of chromosome in which number of chromosome to be kept as the best chromosome among current generation

	int no_generation = 10;		  // number of iteration(main)
	int no_sub_generation = 5;	  // number of iteration(sub, within the unfreeze list group)

	bool use_initial_seed = false;

	bool sga_cyc_subgen = false;

	// --- --- --- --- --- current loop --- --- --- --- --- //
	int current_generation = 0;
	int current_sub_generation = 0;
	bool current_isCycOnlyOptGen = false;

	int Stop_criteria = 5;   // could add the limitation of runtime later according to the paper

	int Elite_size;			// number of chromosome to be kept as the best chromosome among current generation
	int offspring_size;		// number of chromosome to be undergoing GA operation current generation

	// --- --- --- --- --- number of bits in 1 chromosome --- --- --- --- --- //
	int bits_in_chromosome = 8;
	int bits_capability = 255;  // regarding to above bits_in_chromosome 


	// --- --- --- --- --- number of maximum cycle  for VGVC and VGFC --- --- --- --- --- //
	// find out the maximum of cycle time from inital setting, and divided by simulation time, get  number of minmum cycle  
	// defualt in 1, just assume every intersection r FGFC
	int max_no_Cycle = 1;


	//internal use, intersecs_sequence, all translation, generation or random timing plan or the create wholegene wouljd used this list
	std::vector <int> intersecs_sequence;

	///// func could do the input into the following set list 
	//-- - -- - -- - -- - -- - unfreeze_list for freezing the intersecitons while doing SGA -- - -- - -- - -- - -- - //
	std::vector < std::set <int>> unfreeze_list;
	int current_unfreeze_loc = 0;
	//-- - -- - -- - -- - -- - intersections_gp deferentiate out the group for cycle synchronization -- - -- - -- - -- - -- - //
	//key = intersection id itself,   value = intersecs_gps_inf (the cycle would all synchronize with coordinated intersection)
	std::unordered_map<int, intersecs_gps_inf> intersections_gp;

	// the best job from previous generation computation
	I_RunCase* best_ref_opt_job = nullptr;


	// --- --- --- --- --- calculation param after input --- --- --- --- --- //
	// constructor of input_param
	input_param() {
		// the size of the elite to be inherit for the next generation
		this->Elite_size = (int)std::round(this->Populcation_size * this->Elite_rate);
		this->offspring_size = this->Populcation_size - this->Elite_size;
	}
};

}