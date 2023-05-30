#include "Chromosome.h"

#include <iostream>
#include <cassert>
#include <format>

#include "I_RunCase.h"

namespace CTM_GA {

Chromosome::Chromosome(const input_param& sourceData, std::unordered_map<int, intersec_inf>& intersecs_inf)
	:ref_sourceData(sourceData), intersecs_inf(intersecs_inf) {}

Chromosome::Chromosome(const input_param& sourceData, std::unordered_map<int, intersec_inf>& intersecs_inf, std::unordered_map<int, intersec_gene> dec_data_structure)
	:ref_sourceData(sourceData), intersecs_inf(intersecs_inf), dec_data_structure(dec_data_structure) {}

Chromosome::~Chromosome() {}

// get the data inside chromosome object class
const input_param& Chromosome::get_ref_sourceData() const {
	return this->ref_sourceData;
}

std::unordered_map<int, intersec_gene> Chromosome::get_dec_data_structure() const {
	return this->dec_data_structure;
}

void  Chromosome::set_dec_data_structure_by_intersec_id(int intersec_id, intersec_gene&& dec_data_structure) {
	this->dec_data_structure.erase(intersec_id);
	this->dec_data_structure.emplace(intersec_id, std::move(dec_data_structure));
}

void Chromosome::create_initial_seed() {

	for (auto& intersec_id : this->ref_sourceData.intersecs_sequence) {
		// init intersec_gene, for diff intersection
		intersec_gene intersec_gene1{};

		// ```````````````````````````````offset ````````````````````````````````//
		//     only have 1 offset for each intersection
		//--------- intialize the dec_Offset into intersec_gene
		intersec_gene1.dec_Offset = this->intersecs_inf.at(intersec_id).initial_Offset;
		
		// check type
		switch (this->intersecs_inf.at(intersec_id).typeOfplan_per_intersec) {
		case intersec_inf::GAtype_SigCtrl::FGFC:
		{ // start of  FGFC =0
			// ```````````````````````````````Cycle time````````````````````````````````//
			// only have 1 cycle time, int_cycle_vector size = 1
			// False, Cycle time is fixed with original setting
			intersec_gene1.dec_Cycle = { this->intersecs_inf.at(intersec_id).initial_Cycle };

			// ```````````````````````````````phases````````````````````````````````//
			// only have 1 no_cycle , dec_phase_duration size = 1
			// bin2dec convert of phases green
			std::vector<std::vector<int>> phase_duration_collections1;
			for (int no_phase = 0; no_phase < this->intersecs_inf.at(intersec_id).phase_duration_inf_collections.size(); ++no_phase) {
				// No, false, this phase green in this intersection is fixed
				// import directly from from the default value of the certain phase green in phasesinfByintersecs table
				// put phase_gene1 e.g.[X] into phase_duration_collections1 
				phase_duration_collections1.push_back({ this->intersecs_inf.at(intersec_id).phase_duration_inf_collections[no_phase].initial_phase_duration });
			}

			//--------- intialize and add dec_oc_phases_per_intersec to dec_intersections_per_pop
			intersec_gene1.phase_duration_collections = std::move(phase_duration_collections1);

			break;
		}// end of  FGFC =0
		}// end of switch (typeOfplan_intersec_selected)

		set_dec_data_structure_by_intersec_id(intersec_id, std::move(intersec_gene1));
		/// the function actually doing dec_data_structure1.insert(std::pair<int, GA_Module::intersec_gene>(intersec_id, intersec_gene1));
	}
}


// --- --- --- --- --- generate the a population just on the dec_data_structure for the 1st generation--- --- --- --- ---
bool Chromosome::gen_random_timeplan() {
	
	// bool_Check_gene_requirements == false , it would immediately break this func means the generation of this chromsome is failed
	// bool_Check_gene_requirements == true (default) , successfully generated this chromsome
	bool bool_Check_gene_requirements = true;

	// each intersection looping
	for (auto& intersec_id : this->ref_sourceData.intersecs_sequence) {
		//--------- refresh the intersec_gene, for diff intersection
		intersec_gene intersec_gene1{};

		//// ```````````````````````````````check vector of unfreeze_list````````````````````````````````//// 
		bool isUnFreezed = false;
		if (this->ref_sourceData.unfreeze_list.empty() == true) {
			// true, unfreeze_list vector is empty, so all intersection is unfreeze could be changeable
			isUnFreezed = true;
		}
		else if(this->ref_sourceData.current_isCycOnlyOptGen == false){
			// check if this intersection is freezed, 0 means not found , more than or equal to 1 means found intersec_id in that set
			isUnFreezed = this->ref_sourceData.unfreeze_list.at(this->ref_sourceData.current_unfreeze_loc).count(intersec_id);
		}

		//cache current intersec
		auto& currIntersecInf = this->intersecs_inf.at(intersec_id);

		// select the no. of cycle for each crossover run on   0 = 1st cycle, 1= 2nd cycle...
		// by default should be all 0, assume they r all FGFC
		switch (currIntersecInf.typeOfplan_per_intersec) {
			// by default should be all 0, assume they r all FGFC
		case intersec_inf::GAtype_SigCtrl::FGFC:
		{

			//// ```````````````````````````````Cycle time````````````````````````````````////     only have 1 cycle time, int_cycle_vector size = 1
			////............................... check the changeable of Cycle time...............................////	
			int int_cycle = 0;
			std::vector <int > int_cycle_vector{};
			if ((currIntersecInf.Cycle_changeable && isUnFreezed)
				|| (this->ref_sourceData.GA_type == input_param::GAtype::SGA && this->ref_sourceData.current_isCycOnlyOptGen)) {
				//// True, Cycle time is changeable

				// check the intersection group, wheather this intersection id is in the grouo or not
				if (this->ref_sourceData.intersections_gp.count(intersec_id)
					&& intersec_id != this->ref_sourceData.intersections_gp.at(intersec_id).coordinated_intersec_id) {
					// 1, true,intersection group do contain this interseciton id

					// Cycle index in wholegene_ref, it is the coordinated_intersec_id cycle index, no more new string added to wholegene_ref
						// 1 , True, (the cycle would all synchronize with coordinated intersection)
					int_cycle = this->dec_data_structure.at(this->ref_sourceData.intersections_gp.at(intersec_id).coordinated_intersec_id).dec_Cycle[0];
					int_cycle_vector.emplace_back(int_cycle);

				}
				else {
					//first condition false, not in the intersection group or 
					//second condition false, it is the coordinated intersection in which it can be changeable, requires rebuild the location and conversion as well
					// so it could totally ramdomly genrate

					int maxCyc = currIntersecInf.max_Cycle;
					int minCyc = std::max(currIntersecInf.min_Cycle, currIntersecInf.sum_PhaseMinG + currIntersecInf.sum_PhaseInterG);

					//overwrite iff cycle time has changeable interval
					if (this->ref_sourceData.use_initial_seed && currIntersecInf.changableInterval > -1) {
						//take from initial guess
						maxCyc = std::min(maxCyc, (currIntersecInf.initial_Cycle + currIntersecInf.changableInterval));
						minCyc = std::max(minCyc, (currIntersecInf.initial_Cycle - currIntersecInf.changableInterval));
					}

					// range of Max_Cycle_length to   Min_Cycle  (the variation range of the phase that could do the random generation)
					int_cycle = (rand() % (maxCyc - (minCyc - 1))) + minCyc;
					int_cycle_vector.emplace_back(int_cycle);

				}

			}
			else if (this->ref_sourceData.best_ref_opt_job != nullptr) {
				int_cycle_vector.emplace_back(this->ref_sourceData.best_ref_opt_job->get_Chromosome_obj().dec_data_structure.at(intersec_id).dec_Cycle[0]);
			}
			else {
				// this intersec is freezed disallowing variation while isUnFreezed is false
				// False, Cycle time is fixed with original setting
				int_cycle = currIntersecInf.initial_Cycle;
				int_cycle_vector.emplace_back(int_cycle);
			}

			//--------- intialize the dec_Cycle into intersec_gene
			intersec_gene1.dec_Cycle = int_cycle_vector;


			//// ```````````````````````````````offset````````````````````````````````////     only 1 offset accross intersec_gene for each intersection, ignorant on typeOfplan_intersec_selected
			////............................... check the changeable of offsets...............................////
			if (currIntersecInf.Offset_changeable == true && isUnFreezed == true) {
				// this intersec is unfreeze allowing variation while isUnFreezed is true
				
				//--------- intialize the dec_Offset into intersec_gene
				//// True, offset is changeable
				// used previous newly randomly assigned dec_Cycle to limit the variation range of the offset
				intersec_gene1.dec_Offset = rand() % (intersec_gene1.dec_Cycle[0]);

			}
			else if (this->ref_sourceData.best_ref_opt_job != nullptr) {
				intersec_gene1.dec_Offset = (this->ref_sourceData.best_ref_opt_job->get_Chromosome_obj().dec_data_structure.at(intersec_id).dec_Offset);
			}
			else {
				//--------- intialize the dec_Offset into intersec_gene
				//// False, offset is fixed with original setting
				intersec_gene1.dec_Offset = currIntersecInf.initial_Offset;
			}


			//// ```````````````````````````````phases````````````````````````````````////   only have 1 no_cycle , dec_phase_duration size = 1
			//  sum of changeable phases 
			int int_sum_changeable_p = 0;
			// count of changeable phases, for safeguarding 
			int count_changeable_p = 0;
			// sum the min green of that of changeable phases
			int int_sum_changeable_minG = 0;

			// sum the fixed phase either by value or by ratio
			int int_sum_fixed_p = 0;
			// sum the min green of that of fixed phase either by value or by ratio
			int int_sum_fixed_minG = 0;

			//int int_phase;
			int tmp_phase = 0;
			std::vector<int> tmp_phase_per_intersec = {};


			// refresh phase_gen efor each no_phase
			std::vector <int> phase_gene1 = {};
			std::vector<std::vector <int>> phase_duration_collections1 = {};


			for (int no_phase = 0; no_phase < currIntersecInf.phase_duration_inf_collections.size(); ++no_phase) {
				// create that certain amouont of empty phase_gene1 first inside phase_duration_collections1
				phase_duration_collections1.emplace_back(phase_gene1);

				//cache the current phase
				auto& currPhaseInf = currIntersecInf.phase_duration_inf_collections[no_phase];

				int thisMinGreen = currPhaseInf.min_green;
				if (this->ref_sourceData.use_initial_seed && currPhaseInf.changableInterval >= 0) {
					thisMinGreen = std::max(currPhaseInf.min_green,
						(currPhaseInf.initial_phase_duration - currPhaseInf.changableInterval));
				}

				////............................... check the changeable of phase duration...............................////
				// only phase_changeable_by_value and phase_changeable_by_ratio r both in true => could totally randomly generate
				if (currPhaseInf.phase_changeable_by_ratio == true &&
					currPhaseInf.phase_changeable_by_value == true &&
					// this intersec is unfreeze allowing variation while isUnFreezed is true
					isUnFreezed == true) {
					// Yes, True, this phase green in this intersection either changeable by ratio or by value

					// just random number between 0 ~ 100, do the division for getting a random ratio later
					tmp_phase = rand() % 100;
					tmp_phase_per_intersec.emplace_back(tmp_phase);

					// sum up all changeable phase green, for following scaled equation
					int_sum_changeable_p += tmp_phase;
					count_changeable_p++;

					// sum up all changeable minimum green,  for following scaled equation
					int_sum_changeable_minG += thisMinGreen;
				}
				else if (currPhaseInf.phase_changeable_by_ratio == false
					&& currPhaseInf.phase_changeable_by_value == true 
					&& isUnFreezed == true		// this intersec is unfreeze allowing variation while isUnFreezed is true
					|| this->ref_sourceData.current_isCycOnlyOptGen) {					
					// No, false, this phase green in this intersection is fixed by ratio

					if (this->ref_sourceData.current_isCycOnlyOptGen && this->ref_sourceData.best_ref_opt_job != nullptr) {
						//take from best job ref
						tmp_phase = static_cast<int>(std::round(intersec_gene1.dec_Cycle[0] * 
							(this->ref_sourceData.best_ref_opt_job->get_Chromosome_obj().dec_data_structure.at(intersec_id).phase_duration_collections.at(no_phase)[0]
								/ (double)this->ref_sourceData.best_ref_opt_job->get_Chromosome_obj().dec_data_structure.at(intersec_id).dec_Cycle[0])));
					}
					else {
						// the previous newly randomly generated cycle time from the above session
						tmp_phase = static_cast<int>(std::round(intersec_gene1.dec_Cycle[0] * \
							// initial phase duration divided by initial cycle time = fixed ratio
							(currPhaseInf.initial_phase_duration / (double)currIntersecInf.initial_Cycle)));
					}

					// make sure this phase green in this intersection is fixed by ratio, is 100% reach the requirement
					if (tmp_phase > thisMinGreen) {
						// 1, true, tmp_phase > min_green (pass) 
						tmp_phase_per_intersec.emplace_back(tmp_phase);
					}
					else {
						// 0, false, tmp_phase <= min_green (fail) 
						tmp_phase_per_intersec.emplace_back(thisMinGreen);
					}

					// sum up all unchanged phase green, for following scaled equation
					int_sum_fixed_p += tmp_phase;

					// sum up all unchanged minimum green,  for following scaled equation
					int_sum_fixed_minG += thisMinGreen;
				}
				else if (currPhaseInf.phase_changeable_by_ratio == true &&
					currPhaseInf.phase_changeable_by_value == false &&
					// this intersec is unfreeze allowing variation while isUnFreezed is true
					isUnFreezed == true ||
					// or this intersec is freeze disallowing variation while isUnFreezed is false
					isUnFreezed == false) {
					// No, false, this phase green in this intersection is fixed by value

					if (isUnFreezed == false && this->ref_sourceData.best_ref_opt_job != nullptr) {
						//take from best job ref
						tmp_phase = this->ref_sourceData.best_ref_opt_job->get_Chromosome_obj().dec_data_structure.at(intersec_id).phase_duration_collections.at(no_phase)[0];
					}
					else {
						// import directly from from the default value of the certain phase green in phasesinfByintersecs table
						tmp_phase = currPhaseInf.initial_phase_duration;
					}

					tmp_phase_per_intersec.emplace_back(tmp_phase);

					// sum up all unchanged phase green, for following scaled equation
					int_sum_fixed_p += tmp_phase;

					// sum up all unchanged minimum green,  for following scaled equation
					int_sum_fixed_minG += thisMinGreen;
				}
				else
					return false;
			}

			// cal the range of the phase durations is allowed for variation
			int phase_variation_range = intersec_gene1.dec_Cycle[0]
				- (int_sum_changeable_minG + currIntersecInf.sum_PhaseInterG + int_sum_fixed_p);

			// newly assign cycle time -   all intergreen for all phases (sum_PhaseInterG)  - sum of above fixed phase value(including extension and certain min green)
			int unused_range = intersec_gene1.dec_Cycle[0] - (currIntersecInf.sum_PhaseInterG + int_sum_fixed_p);

			std::unordered_map<int, int> t_phaseMin;

			// To scale phases green based on ratio
			for (int no_phase = 0; no_phase < currIntersecInf.phase_duration_inf_collections.size(); ++no_phase) {
				auto& currPhaseInf = currIntersecInf.phase_duration_inf_collections[no_phase];

				////............................... check the changeable of phase duration...............................////
				// only phase_changeable_by_value and phase_changeable_by_ratio r both in true => could totally randomly generate
				if (currPhaseInf.phase_changeable_by_ratio == true &&
					currPhaseInf.phase_changeable_by_value == true &&
					// this intersec is unfreeze allowing variation while isUnFreezed is true
					isUnFreezed == true) {
					// Yes, True, this phase green in this intersection either changeable by ratio or by value

					int thisMinGreen = currPhaseInf.min_green;
					int thisMaxGreen = phase_variation_range + thisMinGreen;
					if (this->ref_sourceData.use_initial_seed && currPhaseInf.changableInterval >= 0) {
						thisMinGreen = std::max(currPhaseInf.min_green,
							(currPhaseInf.initial_phase_duration - currPhaseInf.changableInterval));
					}
					t_phaseMin.emplace(no_phase, thisMinGreen);

					if (int_sum_changeable_p) {
						// make sure all pahses green add up  and their own inter green together, equal to  cycle time
															// range of the phase durations is allowed for variation
						tmp_phase_per_intersec[no_phase] = static_cast<int>(std::round(phase_variation_range
							// the ratio of the changeable (fixed not changeable) phase green
							* (tmp_phase_per_intersec[no_phase] / (double)int_sum_changeable_p))
							// plus their own minium green in order to reach the requirement
							+ thisMinGreen);
					}
					else {
						//ensure non divide by zero
						tmp_phase_per_intersec[no_phase] = static_cast<int>(std::round(phase_variation_range
							// even distribute
							* (1.0 / count_changeable_p))
							// plus their own minium green in order to reach the requirement
							+ thisMinGreen);
					}

					assert(tmp_phase_per_intersec[no_phase] >= 0);
					assert(tmp_phase_per_intersec[no_phase] <= intersec_gene1.dec_Cycle[0]);

					// make sure doesn't exceed changeable Cycle Time
					unused_range -= tmp_phase_per_intersec[no_phase];
				}
				else {
					//// False, this phase green is fixed by either value or ratio
					// nothing do here, coz the true unchaged (fixed not changeable) phase green have been added in above for loop
				}
			}

			//adjust the phases such that no unused left
			if (unused_range != 0) {
				for (auto& minEntry : t_phaseMin) {	//only phases in t_phaseMin can be modified (not fixed by value/ratio)
					//if unused_range is negative, make sure doesn't remove more than minGreen allow
					int canAdd =  std::max((minEntry.second - tmp_phase_per_intersec[minEntry.first]), unused_range);
					tmp_phase_per_intersec[minEntry.first] += canAdd;
					unused_range -= canAdd;
				}
			}
			if (unused_range != 0)
				return false;

			//--------- intialize the dec_phase_duration into phase_duration_collections1 
			for (int no_phase = 0; no_phase < currIntersecInf.phase_duration_inf_collections.size(); ++no_phase)
				phase_duration_collections1[no_phase].emplace_back(tmp_phase_per_intersec[no_phase]);
			intersec_gene1.phase_duration_collections = std::move(phase_duration_collections1);

			// check the crossphase here, 
			// it would only check the unfreeze intersection phases, so its intersection by intersection
			if (isUnFreezed || this->ref_sourceData.current_isCycOnlyOptGen) {
				// this intersec is unfreeze allowing variation while isUnFreezed is true

				bool_Check_gene_requirements = Check_crossPhaseMinGreen(intersec_id, intersec_gene1);

				if (bool_Check_gene_requirements == false) {
					// failed in Check_gene_requirements, so abandon this generation chromosome
					return bool_Check_gene_requirements;
				}
			}

			break;
		} // end of FGFC
		} // end of switch


		// put it tmp data struct into dec_data_structure of its own chromosome object, intersection by intersection
		set_dec_data_structure_by_intersec_id(intersec_id, std::move(intersec_gene1));

	}// end of the intersection id loop

	return bool_Check_gene_requirements;

}



// --- --- --- --- --- create whole gene and the bin location --- --- --- --- ---
void Chromosome::create_wholegene() {
	// clear everything left in the chromosome gene except the dec data
	this->wholegene_ref.clear();
	this->bin_cycle_locations.clear();
	this->bin_offset_locations.clear();
	this->bin_phase_locations.clear();

	//currindex for location
	int currindex = 0;

	// looking at intersection
	// each intersection looping
	for (auto& intersec_id : this->ref_sourceData.intersecs_sequence) {
		auto& currIntersecInf = this->intersecs_inf.at(intersec_id);

		//// ```````````````````````````````check vector of unfreeze_list````````````````````````````````//// 
		bool isUnFreezed = false;
		if (this->ref_sourceData.unfreeze_list.empty() == true) {
			// true, unfreeze_list vector is empty, so all intersection is unfreeze could be changeable
			isUnFreezed = true;
		}
		else if (this->ref_sourceData.current_isCycOnlyOptGen == false) {
			// check if this intersection is freezed
			isUnFreezed = this->ref_sourceData.unfreeze_list.at(this->ref_sourceData.current_unfreeze_loc).count(intersec_id);
		}


		//// ```````````````````````````````offset````````````````````````````````////     only have 1 offset in each intersection
		if (currIntersecInf.Offset_changeable == true && isUnFreezed == true) {
			// make bin value
			this->wholegene_ref.emplace_back(dec2bin(this->dec_data_structure.at(intersec_id).dec_Offset));

			// offset index in wholegene_ref
			this->bin_offset_locations.emplace(intersec_id, currindex);
			currindex++;
		}
		else {
			// false offset is not changeable
			// offset would be the -1 index in wholegene_ref
			this->bin_offset_locations.emplace(intersec_id, -1);
		}


		// select the no. of cycle for each crossover run on   0 = 1st cycle, 1= 2nd cycle...
		// by default should be all 0, assume they r all FGFC
		switch (currIntersecInf.typeOfplan_per_intersec) {
			// by default should be all 0, assume they r all FGFC
		case intersec_inf::GAtype_SigCtrl::FGFC:
		{ // start of  FGFC =0

			std::vector<int> t_cycle_loc;
			//// ```````````````````````````````Cycle time````````````````````````````````////     only have 1 cycle time, int_cycle_vector size = 1
			if ((currIntersecInf.Cycle_changeable && isUnFreezed)
				|| (this->ref_sourceData.GA_type == input_param::GAtype::SGA && this->ref_sourceData.current_isCycOnlyOptGen)) {
				//// True, Cycle time is changeable

				// check the intersection group, wheather this intersection id is in the grouo or not
				if (this->ref_sourceData.intersections_gp.count(intersec_id)
					&& intersec_id != this->ref_sourceData.intersections_gp.at(intersec_id).coordinated_intersec_id) {
					// 1, true,intersection group do contain this interseciton id
					// 1 , True, (the cycle would all synchronize with coordinated intersection)
					// 2 false, it is the coordinated intersection in which it can be changeable, requires rebuild the location and conversion as well

					// Cycle index in wholegene_ref, it is the coordinated_intersec_id cycle index, no more new string added to wholegene_ref
					t_cycle_loc.emplace_back(this->bin_cycle_locations.at(this->ref_sourceData.intersections_gp.at(intersec_id).coordinated_intersec_id)[0]);
					// so no add one of the currindex
				}
				else {
					// 0, flase, current intersection is not included in the intersection group

					int maxCyc = currIntersecInf.max_Cycle;
					int minCyc = std::max(currIntersecInf.min_Cycle, currIntersecInf.sum_PhaseMinG + currIntersecInf.sum_PhaseInterG);

					//overwrite iff cycle time has changeable interval
					if (this->ref_sourceData.use_initial_seed && currIntersecInf.changableInterval > -1) {
						//take from initial guess
						maxCyc = std::min(maxCyc, (currIntersecInf.initial_Cycle + currIntersecInf.changableInterval));
						minCyc = std::max(minCyc, (currIntersecInf.initial_Cycle - currIntersecInf.changableInterval));
					}

					// conversion from the decimal value to string by the original scale equation
					// so the string results acutally keep the genuine ratio features
					this->wholegene_ref.emplace_back(dec2bin((int)std::round(
						((double)this->dec_data_structure.at(intersec_id).dec_Cycle[0] - minCyc)
						* (this->ref_sourceData.bits_capability)
						/ ((double)maxCyc - minCyc)
					)));

					// Cycle index in wholegene_ref for the coordinated intersection
					t_cycle_loc.emplace_back(currindex);
					currindex++;
				}
			}
			else {
				// this intersec is freezed disallowing variation while isUnFreezed is false
				// false Cycle is not changeable
				// Cycle would always be the 1 index in wholegene_ref
				t_cycle_loc.emplace_back(-1);
			}
			// put tmp vector of loc of cycle to bin_cycle_locations
			this->bin_cycle_locations.emplace(intersec_id, std::move(t_cycle_loc));


			std::vector<std::vector<int>> phases_loc;   // for multiple phase 
			std::vector<int> phase_loc;   // for 1 phase 
			//// ```````````````````````````````phases````````````````````````````````////   only have 1 no_cycle , dec_phase_duration size = 1
			for (auto phase_no = 0; phase_no < currIntersecInf.phase_duration_inf_collections.size(); phase_no++) {
				phase_loc.clear();

				// only phase_changeable_by_value and phase_changeable_by_ratio r both in true => could totally randomly generate
				if (currIntersecInf.phase_duration_inf_collections[phase_no].phase_changeable_by_ratio == true &&
					currIntersecInf.phase_duration_inf_collections[phase_no].phase_changeable_by_value == true &&
					// this intersec is unfreeze allowing variation while isUnFreezed is true
					isUnFreezed == true) {
					// Yes, True, this phase green in this intersection either changeable by ratio or by value

					this->wholegene_ref.emplace_back(dec2bin(this->dec_data_structure.at(intersec_id).phase_duration_collections[phase_no][0]));

					// phase index in wholegene_ref
					phase_loc.emplace_back(currindex);
					currindex++;
				}
				else {
					// this intersec is freezed disallowing variation while isUnFreezed is false

					// false phase is not changeable
					// phase would always be the 1 index in wholegene_ref
					phase_loc.emplace_back(-1);
				}
				phases_loc.emplace_back(phase_loc);

			}
			// put tmp vector of loc of cycle to bin_cycle_locations
			this->bin_phase_locations.emplace(intersec_id, std::move(phases_loc));
			break;
		}// end of  FGFC =0

		}// end of switch (typeOfplan_intersec_selected)
	}// end of the intersection id loop
}

// --- --- --- --- --- generate set of crossover point locations for one crossover and mutation on a pair of chromosome--- --- --- --- ---
// usde for the crossoverBypt_mutationByrate() func
std::set<int> Chromosome::generate_rand_crossover_point() {
	int num_ele = this->wholegene_ref.size();
	if (static_cast<int>(num_ele * this->ref_sourceData.crossover_multiplier) >
			(num_ele * this->ref_sourceData.bits_in_chromosome)) {
		throw std::runtime_error("too many crossover point is exceeding the limit");
	}

	// gen non-repeating crossover pts
	std::set<int> out;
	for (int pt_no = 0; pt_no < std::round(num_ele * this->ref_sourceData.crossover_multiplier); pt_no++) {
		int tmp_pt;
		do {
			tmp_pt = rand() % (num_ele * this->ref_sourceData.bits_in_chromosome);
		} while (out.count(tmp_pt));
		out.insert(tmp_pt);
	}
	return out;
}

// --- --- --- --- --- according to bin location --- --- --- --- ---
bool Chromosome::wholegene_operation2dec_data() {

	bool bool_Check_gene_requirements = true;

	// looking at intersection
	// each intersection looping

	for (auto& intersec_id : this->ref_sourceData.intersecs_sequence) {
		auto& currIntersecInf = this->intersecs_inf.at(intersec_id);

		//// ```````````````````````````````check vector of unfreeze_list````````````````````````````````//// 
		bool isUnFreezed = false;
		if (this->ref_sourceData.unfreeze_list.empty() == true) {
			// true, unfreeze_list vector is empty, so all intersection is unfreeze could be changeable
			isUnFreezed = true;
		}
		else if (this->ref_sourceData.current_isCycOnlyOptGen == false) {
			// check if this intersection is freezed, 0 means not found , more than or equal to 1 means found intersec_id in that set
			isUnFreezed = this->ref_sourceData.unfreeze_list.at(this->ref_sourceData.current_unfreeze_loc).count(intersec_id);
		}

		// select the no. of cycle for each crossover run on   0 = 1st cycle, 1= 2nd cycle...
		// by default should be all 0, assume they r all FGFC
		switch (currIntersecInf.typeOfplan_per_intersec) {
			// by default should be all 0, assume they r all FGFC
		case intersec_inf::GAtype_SigCtrl::FGFC:
		{ // start of  FGFC =0

			////............................... check the changeable of Cycle time...............................////

			//// True, Cycle time is changeable
			if ((currIntersecInf.Cycle_changeable && isUnFreezed)
				|| (this->ref_sourceData.GA_type == input_param::GAtype::SGA && this->ref_sourceData.current_isCycOnlyOptGen)) {

				// tmp for scaled cycle from binary string to decimal integer
				std::vector <int> scaled_int_cycle;

				if (this->ref_sourceData.intersections_gp.count(intersec_id)
					&& intersec_id != this->ref_sourceData.intersections_gp.at(intersec_id).coordinated_intersec_id) {

					// true, intersection is in the group but the current intersection is not the coordinated intersection
					// directly copy the exact decimal value of cycle time from the coordinated intersection
					scaled_int_cycle.emplace_back(this->dec_data_structure.at(this->ref_sourceData.intersections_gp.at(intersec_id).coordinated_intersec_id).dec_Cycle[0]);

					// replace the updated real value back to string bin_oc_phases_per_intersec
					this->dec_data_structure.at(intersec_id).dec_Cycle[0] = scaled_int_cycle[0];

				} // end of intersections_gp check for ture intersection contains in the group
				else {
					// 0 , false, intersection group dosent contain this interseciton id\
					// 1 , True, (the cycle would all synchronize with coordinated intersection)
					// cycle time could change totally randomly but those r in the same group would have the same cycle time, due to minimum_maxCycle, maximum_minCycle, and that of wholegene_ref r the same
					// bin_cycle_locations point to the same loaction in the wholegene_operation

					int cycle_min = std::max(currIntersecInf.min_Cycle, currIntersecInf.sum_PhaseMinG + currIntersecInf.sum_PhaseInterG);
					int cycle_max = currIntersecInf.max_Cycle;  // Max_Cycle_length

					//overwrite iff cycle time has changeable interval
					if (this->ref_sourceData.use_initial_seed && currIntersecInf.changableInterval > -1) {
						//take from initial guess
						cycle_max = std::min(cycle_max, (currIntersecInf.initial_Cycle + currIntersecInf.changableInterval));
						cycle_min = std::max(cycle_min, (currIntersecInf.initial_Cycle - currIntersecInf.changableInterval));
					}

					// cycle time could change totally randomly
					scaled_int_cycle.emplace_back(static_cast<int> (std::round(bin2dec(this->wholegene_ref.at(this->bin_cycle_locations.at(intersec_id)[0]))
						* (((double)cycle_max - cycle_min) / static_cast<double>(this->ref_sourceData.bits_capability)) + cycle_min)));
					// replace the updated real value back to string bin_oc_phases_per_intersec
					this->dec_data_structure.at(intersec_id).dec_Cycle[0] = scaled_int_cycle[0];

				} // end of intersections_gp check for false intersection not in the group
			}// end of Cycle_changeable and isUnFreezed check for true



			//// ```````````````````````````````offset````````````````````````````````////     only have 1 offset in each intersection
			// no need to check the intersection group, offset is irrelated to the intersection group, could be totaly random, as long as it not exceeding the cycle time
			if (this->intersecs_inf.at(intersec_id).Offset_changeable == true &&
				// this intersec is unfreeze allowing variation while isUnFreezed is true
				isUnFreezed == true) {
				//// True, offset is changeable

				// copy out a session of string from the wholegene_operation
				std::string tmp_bin_offset;
				tmp_bin_offset = this->wholegene_ref.at(this->bin_offset_locations.at(intersec_id));

				// direct conversion from the string to the value
				this->dec_data_structure.at(intersec_id).dec_Offset = bin2dec(tmp_bin_offset) % this->dec_data_structure.at(intersec_id).dec_Cycle[0];

			}// end of Offset_changeable and isUnFreezed check for true
			else {
				//// False, Offset wouldnt undergoes crossover and mutation, its not changeable either becuase of intersection is freezed or changeable bool is false
				this->dec_data_structure.at(intersec_id).dec_Offset = this->dec_data_structure.at(intersec_id).dec_Offset % this->dec_data_structure.at(intersec_id).dec_Cycle[0];

			}// end of Offset_changeable and isUnFreezed check for False



			//// ```````````````````````````````phases````````````````````````````````////   only have 1 no_cycle , dec_phase_duration size = 1
			// phases wouljd undergo readjustment no matter what, so cycle_has_changes is not used in this part in order to avoid error  to make sure every phase do reach requirement,
			// but only afraid of the randomly generate cycle time is too small, which makes some phase become negative value

			//  sum of changeable phases 
			int int_sum_changeable_p = 0;
			// count of changeable phases, for safeguarding 
			int count_changeable_p = 0;
			// sum the min green of that of changeable phases
			int int_sum_changeable_minG = 0;

			// sum the fixed phase either by value or by ratio
			int int_sum_fixed_p = 0;
			// sum the min green of that of fixed phase either by value or by ratio
			int int_sum_fixed_minG = 0;

			//int int_phase;
			int tmp_phase = 0;
			std::vector<int> tmp_phase_per_intersec{};


			// refresh phase_gen efor each no_phase
			std::vector <int> phase_gene1{};
			std::vector<std::vector <int>> phase_duration_collections1{};

			for (int no_phase = 0; no_phase < currIntersecInf.phase_duration_inf_collections.size(); ++no_phase) {
				// create that certain amouont of empty phase_gene1 first inside phase_duration_collections1
				phase_duration_collections1.emplace_back(phase_gene1);

				//cache the current phase
				auto& currPhaseInf = currIntersecInf.phase_duration_inf_collections[no_phase];

				int thisMinGreen = currPhaseInf.min_green;
				if (this->ref_sourceData.use_initial_seed && currPhaseInf.changableInterval >= 0) {
					thisMinGreen = std::max(currPhaseInf.min_green,
						(currPhaseInf.initial_phase_duration - currPhaseInf.changableInterval));
				}

				////............................... check the changeable of phase duration...............................////
				// only phase_changeable_by_value and phase_changeable_by_ratio r both in true => could totally randomly generate
				if (currPhaseInf.phase_changeable_by_ratio == true &&
					currPhaseInf.phase_changeable_by_value == true &&
					// this intersec is unfreeze allowing variation while isUnFreezed is true
					isUnFreezed == true) {
					// this intersec is unfreezed allowing variation while isUnFreezed is true
					// Yes, True, this phase green in this intersection either changeable by ratio or by value

					// getting the string offset index starting point
					int phase_index;
					phase_index = this->bin_phase_locations.at(intersec_id).at(no_phase)[0];

					// copy out a session of string from the wholegene_operation
					std::string tmp_bin_phase;
					tmp_bin_phase = this->wholegene_ref.at(phase_index);

					//  translate binary offset into decimal integer
					int tmp_phase = bin2dec(tmp_bin_phase);

					tmp_phase_per_intersec.emplace_back(tmp_phase);

					// sum up all changeable phase green, for following scaled equation
					int_sum_changeable_p += tmp_phase;
					++count_changeable_p;

					// sum up all changeable minimum green,  for following scaled equation
					int_sum_changeable_minG += thisMinGreen;
				}
				else if (currPhaseInf.phase_changeable_by_ratio == false &&
					currPhaseInf.phase_changeable_by_value == true &&
					// this intersec is unfreeze allowing variation while isUnFreezed is true
					isUnFreezed == true
					|| this->ref_sourceData.current_isCycOnlyOptGen) {
					// this intersec is freezed disallowing variation while isUnFreezed is false
					// or isUnFreezed is true, phase_changeable_by_ratio == false , this phase green in this intersection is fixed by ratio

									// the previous newly randomly generated cycle time from the above session
					tmp_phase = static_cast<int>(std::round(this->dec_data_structure.at(intersec_id).dec_Cycle[0] * \
						// initial phase duration divided by initial cycle time = fixed ratio
						(currPhaseInf.initial_phase_duration / (double)currIntersecInf.initial_Cycle)));

					// make sure this phase green in this intersection is fixed by ratio, is 100% reach the requirement
					if (tmp_phase > thisMinGreen) {
						// 1, true, tmp_phase > min_green (pass) 
						tmp_phase_per_intersec.emplace_back(tmp_phase);
					}
					else {
						// 0, false, tmp_phase <= min_green (fail) 
						tmp_phase_per_intersec.emplace_back(thisMinGreen);
					}

					// sum up all unchanged phase green, for following scaled equation
					int_sum_fixed_p += tmp_phase;

					// sum up all unchanged minimum green,  for following scaled equation
					int_sum_fixed_minG += thisMinGreen;
				}
				///////  ------   coz this is fixed as initial setting already no need to reassign into the dec_data_structure again
				else if (currPhaseInf.phase_changeable_by_ratio == true &&
					currPhaseInf.phase_changeable_by_value == false &&
					// this intersec is unfreeze allowing variation while isUnFreezed is true
					isUnFreezed == true ||
					// or this intersec is freeze disallowing variation while isUnFreezed is false
					isUnFreezed == false) {
					// this intersec is freezed disallowing variation while isUnFreezed is false
					// or isUnFreezed is true, phase_changeable_by_value == false , this phase green in this intersection is fixed by value

									// import directly from from the default value of the certain phase green in phasesinfByintersecs table
					tmp_phase = currPhaseInf.initial_phase_duration;
					tmp_phase_per_intersec.emplace_back(tmp_phase);

					// sum up all unchanged phase green, for following scaled equation
					int_sum_fixed_p += tmp_phase;

					// sum up all unchanged minimum green,  for following scaled equation
					int_sum_fixed_minG = thisMinGreen;
				}
				else
					return false;
			}


			// cal the range of the phase durations is allowed for variation
			int phase_variation_range = this->dec_data_structure.at(intersec_id).dec_Cycle[0]
				- (currIntersecInf.sum_PhaseInterG + int_sum_changeable_minG + int_sum_fixed_p);

			// newly assign cycle time -   all intergreen for all phases (sum_PhaseInterG)  - sum of above fixed phase value(including extension and certain min green)
			int unused_range = this->dec_data_structure.at(intersec_id).dec_Cycle[0] - (currIntersecInf.sum_PhaseInterG + int_sum_fixed_p);

			std::unordered_map<int, int> t_phaseMin;

			// check of unused_range having enough room
			if (unused_range < int_sum_changeable_minG) {
				// unused_range rest of the green variation range is not enough for distruibuting to the changeable phases
				// therefore forced to stop the func of generation chromosome
				bool_Check_gene_requirements = false;
				return bool_Check_gene_requirements;
			}

			// To scale phases green based on ratio
			for (int no_phase = 0; no_phase < currIntersecInf.phase_duration_inf_collections.size(); ++no_phase) {
				auto& currPhaseInf = currIntersecInf.phase_duration_inf_collections[no_phase];

				////............................... check the changeable of phase duration...............................////
				// only phase_changeable_by_value and phase_changeable_by_ratio r both in true => could totally randomly generate
				if (currPhaseInf.phase_changeable_by_ratio == true &&
					currPhaseInf.phase_changeable_by_value == true &&
					// this intersec is unfreeze allowing variation while isUnFreezed is true
					isUnFreezed == true) {
					// Yes, True, this phase green in this intersection either changeable by ratio or by value

					int thisMinGreen = currPhaseInf.min_green;
					int thisMaxGreen = phase_variation_range;
					if (this->ref_sourceData.use_initial_seed && currPhaseInf.changableInterval >= 0) {
						thisMinGreen = std::max(currPhaseInf.min_green,
							(currPhaseInf.initial_phase_duration - currPhaseInf.changableInterval));
					}

					if (int_sum_changeable_p) {
						// make sure all pahses green add up  and their own inter green together, equal to  cycle time
															// range of the phase durations is allowed for variation
						tmp_phase_per_intersec[no_phase] = static_cast<int>(std::round(phase_variation_range
							// the ratio of the changeable (fixed not changeable) phase green
							* (tmp_phase_per_intersec[no_phase] / (double)int_sum_changeable_p))
							// plus their own minium green in order to reach the requirement
							+ thisMinGreen);
					}
					else {
						//ensure non divide by zero
						tmp_phase_per_intersec[no_phase] = static_cast<int>(std::round(phase_variation_range
							// even distribute
							* (1.0 / count_changeable_p))
							// plus their own minium green in order to reach the requirement
							+ thisMinGreen);
					}

					assert(tmp_phase_per_intersec[no_phase] >= 0);
					assert(tmp_phase_per_intersec[no_phase] <= this->dec_data_structure.at(intersec_id).dec_Cycle[0]);

					// make sure doesn't exceed changeable Cycle Time
					unused_range -= tmp_phase_per_intersec[no_phase];
				}
				else {
					//// False, this phase green is fixed by either value or ratio
					// nothing do here, coz the true unchaged (fixed not changeable) phase green have been added in above for loop
				}
			}

			//adjust the phases such that no unused left
			if (unused_range != 0) {
				for (auto& minEntry : t_phaseMin) {	//only phases in t_phaseMin can be modified (not fixed by value/ratio)
					//if unused_range is negative, make sure doesn't remove more than minGreen allow
					int canAdd = std::max((minEntry.second - tmp_phase_per_intersec[minEntry.first]), unused_range);
					tmp_phase_per_intersec[minEntry.first] += canAdd;
					unused_range -= canAdd;
				}
			}
			if (unused_range != 0)
				return false;

			//--------- intialize the dec_phase_duration into phase_duration_collections1 
			for (int no_phase = 0; no_phase < currIntersecInf.phase_duration_inf_collections.size(); ++no_phase)
				phase_duration_collections1[no_phase].emplace_back(tmp_phase_per_intersec[no_phase]);

			// clear for the next round diff intersection
			tmp_phase_per_intersec.clear();

			//--------- replace the new phase_duration_collections into the dec_data_structure phase part
			this->dec_data_structure.at(intersec_id).phase_duration_collections.clear();
			this->dec_data_structure.at(intersec_id).phase_duration_collections = phase_duration_collections1;

			// check the crossphase here, 
			// it would only check the unfreeze intersection phases, so its intersection by intersection
			if (isUnFreezed == true) {
				// this intersec is unfreeze allowing variation while isUnFreezed is true
				bool_Check_gene_requirements = Check_crossPhaseMinGreen(intersec_id, this->dec_data_structure.at(intersec_id));

				if (bool_Check_gene_requirements == false) {
					// this equation would jump into the end, coz bool_Check_gene_requirements r false (fail in generating this chromosome)
					return bool_Check_gene_requirements;
				}
			}

			break;
		}// end of  FGFC =0
		}// end of switch (typeOfplan_intersec_selected)

	}// end of the intersection id loop

	return bool_Check_gene_requirements;
}



// --- --- --- --- --- GA processing --- --- --- --- ---
	// selected the chromosome already the input to here
////............................... crossover by rate and mutation by rate ...............................////
// do crossover by rate and mutation by rate on two selected chromosome in the  jobs_ptr_collections pool
// it could be repeatedly called if the chromosome is failed in check validity
void Chromosome::crossoverByrate_mutationByrate(Chromosome& c) {
	bool swapFlag = false;  // true = do swap    // flase = no swap

	// looping selement wise inside in the vecotr of string of wholegene
	for (size_t str_index = 0; str_index < this->wholegene_ref.size(); str_index++) {
		// looping character wise inside each element in the vecotr of string of wholegene
		for (int ch_loc = 0; ch_loc < this->wholegene_ref[str_index].size(); ch_loc++) {

			////............................... crossover by rate ...............................////
			// random random generate a double float 
			if ((rand() % 100) < this->ref_sourceData.Crossover_rate * 100)
				swapFlag = !swapFlag;
			//swap iff swap flag is active
			if (swapFlag) {
				char temp = this->wholegene_ref[str_index][ch_loc];
				this->wholegene_ref[str_index][ch_loc] = c.wholegene_ref[str_index][ch_loc];
				c.wholegene_ref[str_index][ch_loc] = temp;
			}

			////............................... mutation by rate ...............................////
			// random generate a percentage if its smaller than mutation rate then mutation occurs on wholegene_ref chromosome1
			if ((rand() % 100) < this->ref_sourceData.Mutation_rate * 100)
				bin_flip(this->wholegene_ref[str_index][ch_loc]);
			if ((rand() % 100) < this->ref_sourceData.Mutation_rate * 100)
				bin_flip(c.wholegene_ref[str_index][ch_loc]);
		}// end of ch location index loop
	}// end of string index loop
}

void Chromosome::crossoverBypt_mutationByrate(Chromosome& c) {
	////............................... crossover by crossover point ...............................////
	// generate set of crossover point locations for crossover position on a pair of chromosome
	auto crossover_pts = generate_rand_crossover_point();
	int currPt = 0;			// keeps track if it is a crossover pt
	bool swapFlag = false;  // true = do swap    // flase = no swap

	// looping selement wise inside in the vecotr of string of wholegene
	for (size_t str_index = 0; str_index < this->wholegene_ref.size(); str_index++) {
		// looping character wise inside each element in the vecotr of string of wholegene
		for (int ch_loc = 0; ch_loc < this->wholegene_ref[str_index].size(); ch_loc++) {
			//flip flag if it is crossing point
			if (crossover_pts.count(currPt))
				swapFlag = !swapFlag;
			//swap iff swap flag is active
			if (swapFlag) {
				char temp = this->wholegene_ref[str_index][ch_loc];
				this->wholegene_ref[str_index][ch_loc] = c.wholegene_ref[str_index][ch_loc];
				c.wholegene_ref[str_index][ch_loc] = temp;
			}

			////............................... mutation by rate ...............................////
			// random generate a percentage if its smaller than mutation rate then mutation occurs on wholegene_ref chromosome1
			if ((rand() % 100) < this->ref_sourceData.Mutation_rate * 100)
				bin_flip(this->wholegene_ref[str_index][ch_loc]);
			if ((rand() % 100) < this->ref_sourceData.Mutation_rate * 100)
				bin_flip(c.wholegene_ref[str_index][ch_loc]);

			// update the character position among the wholegene
			++currPt;
		}// end of ch location index loop
	}// end of string index loop
}

// for a single letter character to flip from 1 to 0, 0 to 1
// this funciton is inside the mutation function
char Chromosome::bin_flip(char& ch) {
	ch = (ch == '1') ? '0' : '1';
	return ch;
}

int Chromosome::bin2dec(std::string bin) {
	int out = 0;
	auto itr = bin.rbegin();
	for (int i = 0; i < bin.size(); i++) {
		out += (int)(pow(2, i) * ((double)*itr - '0'));
		itr++;
	}
	return out;
}

std::string Chromosome::dec2bin(int dec) {
	std::string out(this->ref_sourceData.bits_in_chromosome, '0');
	for (int i = this->ref_sourceData.bits_in_chromosome - 1; i >= 0; i--) {
		if (dec % 2)
			++(out.at(i));
		dec /= 2;
	}
	return out;
}

bool Chromosome::Check_crossPhaseMinGreen(const int& intersec_id, const intersec_gene& intersec_obj) {

	// firstly assume all phases in this intersection do reach the crossphase requirement whether it has the crossphase or not
	bool pass = true;

	// --- --- --- --- --- --- --- --- check of  PhaseidList_CrossPhasesMinGreen std::pair first : list of phases index, second : crossphase min green  --- --- --- --- --- --- --- ---  //
	if (this->intersecs_inf.at(intersec_id).PhaseidList_CrossPhasesMinGreen.empty()) {
		// this intersection have no crossphase with empty vector here
		// do nothing
	}
	else {
		// this phase duration in this intersection have crossphase

		// define the temp value for following crosspahse value checking 
		int current_check_phase = 0, previous_phase_durations, previous_intergreens;
		int crossphase_druation;

		for (std::size_t no_cycle = 0; no_cycle < intersec_obj.dec_Cycle.size(); ++no_cycle) {
			previous_phase_durations = 0;
			previous_intergreens = 0;

			//  std::pair first : list of phases index, second : crossphase min green  =====> this is for the crossphase only, the vector size depends on how many crossphase u have
			for (auto Crossphase1 : this->intersecs_inf.at(intersec_id).PhaseidList_CrossPhasesMinGreen) {

				// list of crossphases index (eg. <2,3,4> or <4,0>), looping each crossphase index
				for (auto crossphase_index_it = Crossphase1.first.begin(); crossphase_index_it != Crossphase1.first.end(); ++crossphase_index_it) {

					if (*crossphase_index_it == Crossphase1.first.back()) {
						// crossphase_index loop to the last element in which the real crossphase value checking starts
						// if loop to the last element phase of the crossphase, dont need to add up,saved that of phase duration into current_check_phase
						current_check_phase = intersec_obj.phase_duration_collections[*crossphase_index_it][no_cycle];

					}
					else {
						//list of crossphases index (eg. <2,3,4> or <4,0>) the crossphase_index was still havent reached the last element

						//  add up previous phase of the crossphase
						previous_phase_durations += intersec_obj.phase_duration_collections[*crossphase_index_it][no_cycle];
						//  add up previous inter green of the crossphase
						previous_intergreens += this->intersecs_inf.at(intersec_id).phase_duration_inf_collections[*crossphase_index_it].intergreen;

					}

				}

				crossphase_druation = current_check_phase + previous_phase_durations + previous_intergreens;

				if (crossphase_druation < Crossphase1.second) {
					// true, pass to here means crossPhaseMinGreen is failed for this no_cycle or for this intersection; failed
					pass = false;
					return pass;
				}
			}
		}
	}
	return pass;
}

void Chromosome::print_dec_data_structure() {
	for (auto intersec_gene : this->dec_data_structure) {
		std::cout << "\n intersection_id: " << intersec_gene.first << "\t";
		std::cout << "Cycle: " << intersec_gene.second.dec_Cycle[0] << "\t";
		std::cout << "Offset: " << intersec_gene.second.dec_Offset << "\t";

		std::cout << "phase : [";
		for (auto& phase : intersec_gene.second.phase_duration_collections)
			std::cout << phase[0] << "\t";
		std::cout << " ]";
	}
	std::cout << "\n" << std::endl; ;
}

}	//namespace CTM_GA 