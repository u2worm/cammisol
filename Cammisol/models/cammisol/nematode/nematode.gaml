/**
* Name: nematode
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model nematode

import "../microbes/microbes.gaml"

global {
	float nematode_CO2_emissions <- 0.0#gram;
	float nematode_predation_rate <- 1e5 * 1e-12 #gram/#day; // About 100k bacterias (1e-12#gram) each day
	float nematode_C_N <- 6.0;
	float nematode_C_P <- 47.0;
	float nematode_CUE <- 0.16;
	float nematode_sleep_threshold <- 0.1;
	float nematode_rejected_microbe_necromass_rate <- 0.5;
	// [Initialisation] : taille (côté) minimale d'un pore accessible aux nématodes,
	// pour le placement initial ET les déplacements. 0 = aucune exclusion.
	float taille_min_pore_nematode <- 0.0 #m;
}

species Nematode schedules:[]
{
	// quantity of C N P
//	float C <- (470#gram * 10^-6); // TODO: source?
//	float N <- (470#gram * 10^-6)/C_N;
//	float P <- (470#gram * 10^-6)/C_P;
	
	float stomack_C;
	float stomack_N;
	float stomack_P;
	
	float requested_C;
	float perceived_C;
	
	PoreParticle current_pore;
	
	bool awake <- true;
//	float predation_rate <- 100*0.5e-9#gram/#mn;
//	float predation_rate <- 0.00001#gram/#day; // 0.5E-9#gram/#mn;
	
	//un nematode mange entre 10K et 100k bacteries par jour -> moins de 10k il dort (ref article de colman)
	
	action life() {
	    requested_C <- nematode_predation_rate * local_step;
	    
	    do move();
	    perceived_C <- perceive(); // apres avoir bougé
	    
	    awake <- requested_C > 0 and (perceived_C / requested_C) >= nematode_sleep_threshold;
	    if awake { do predation(); do respirate(); do anabolize(); }
	}
	
	float perceive(){
		return sum(current_pore.populations collect(each.C + each.cytosol_C));
	}
	
	action move()
	{
		// : un nématode n'entre pas dans un pore plus petit que taille_min_pore_nematode
		list<PoreParticle> neighbors_pores <- current_pore.pore_neighbors where (each.size >= taille_min_pore_nematode);
		
		//write "pore courant" + current_pore + ": x: " + current_pore.grid_x + ", y: " + current_pore.grid_y + ", z:" + current_pore.grid_z;
		
		//write "pores voisins : " + length(neighbors_pores) + " \n ";
		//loop n over: neighbors_pores{
		//	write "  voisin " + n + " : x: " + n.grid_x + ", y: " + n.grid_y + ", z: " + n.grid_z;
		//}
		
		float max_perceive <- perceived_C;
		PoreParticle most_attractive <- current_pore;
		ask neighbors_pores
		{
			float current_perceived_C <- sum(populations collect (each.C + each.cytosol_C));
			if(current_perceived_C > max_perceive) {
				max_perceive <- current_perceived_C;
				most_attractive <- self;
			}
		}
		current_pore <- most_attractive;
		location <- any_location_in(current_pore); //TODO a changer
		
		//write "nouveau pore" + current_pore + "; x: " + current_pore.grid_x + ", y: " + current_pore.grid_y + ", z:" + current_pore.grid_z;
	}
	
	action predation()
 	{
		float total_C_predated <- min([perceived_C, requested_C]);
		//write "total_C_predated : " + total_C_predated + "(" + perceived_C +"," + requested_C + ")" ;
		// Total C predated is computed from global perceptions, but N and P predated must be computed from the local state of each population.
		float total_N_predated <- 0.0;
		float total_P_predated <- 0.0;
		
		ask current_pore.populations
		{	
			// perceived_C = sum(C + cytocol_C) on all the populations. So, if total_C_predated = perceived_C,
			// all the C is consumed from each population. Else if total_C_predated = requested_C, an equal proportion
			// of total_C_predated/perceived_C is consumed from each population so that the total C consumed is equal to
			// requested_C.
			float total_C_predated_from_this_population <- (total_C_predated / myself.perceived_C) * (C + cytosol_C); //pourquoi cette formule ?
			
			//write "total_C_predated_from_this_population " + total_C_predated_from_this_population; 
			
			if (total_C_predated_from_this_population > 0.0) {
				float C_predated_from_structural <- C/(C+cytosol_C) * total_C_predated_from_this_population;
				float N_predated_from_structural <- C_predated_from_structural * N/C;
				float P_predated_from_structural <- C_predated_from_structural * P/C;
				
				//write "C_predated_from_structural" + C_predated_from_structural;
				//write "N_predated_from_structural" + N_predated_from_structural;
				//write "P_predated_from_structural" + P_predated_from_structural;
							
				C <- C - C_predated_from_structural;
				N <- N - N_predated_from_structural;
				P <- P - P_predated_from_structural;
				
				total_N_predated <- total_N_predated + N_predated_from_structural;
				total_P_predated <- total_P_predated + P_predated_from_structural;
				
				float C_predated_from_cytosol <- total_C_predated_from_this_population - C_predated_from_structural;
				
				//write "C_predated_from_cytosol" + C_predated_from_cytosol;
				
				if(C_predated_from_cytosol > 0.0) {
					float N_predated_from_cytosol <- C_predated_from_cytosol * cytosol_N / cytosol_C;	
					float P_predated_from_cytosol <- C_predated_from_cytosol * cytosol_P / cytosol_C;
					
					//write "N_predated_from_cytosol" + N_predated_from_cytosol;
					//write "P_predated_from_cytosol" + P_predated_from_cytosol;
						
					cytosol_C <- cytosol_C - C_predated_from_cytosol;
					cytosol_N <- cytosol_N - N_predated_from_cytosol;
					cytosol_P <- cytosol_P - P_predated_from_cytosol;
					
					total_N_predated <- total_N_predated + N_predated_from_cytosol;
					total_P_predated <- total_P_predated + P_predated_from_cytosol;
				}
			}
		}
		
		float microbe_necromass_C <- total_C_predated*nematode_rejected_microbe_necromass_rate;
		float microbe_necromass_N <- total_N_predated*nematode_rejected_microbe_necromass_rate;
		float microbe_necromass_P <- total_P_predated*nematode_rejected_microbe_necromass_rate;
		
		do reject_to_dom(microbe_necromass_C, microbe_necromass_N, microbe_necromass_P); //Rejeté dans le dam (changement)
		
		stomack_C <- stomack_C + total_C_predated - microbe_necromass_C;
		stomack_N <- stomack_N + total_N_predated - microbe_necromass_N;
		stomack_P <- stomack_P + total_P_predated - microbe_necromass_P;
 	}
	
	action respirate()
	{
		float C_respiration <- (1-nematode_CUE) * stomack_C;
		stomack_C <- stomack_C  - C_respiration;
		nematode_CO2_emissions <- nematode_CO2_emissions + C_respiration;
	}
	
	action anabolize()
	{
		float C_assimilated <- min(stomack_C, stomack_N * nematode_C_N, stomack_P * nematode_C_P);
		stomack_C <- stomack_C  - C_assimilated; 
		
		float N_assimilated <- C_assimilated/nematode_C_N;
		float P_assimilated <- C_assimilated/nematode_C_P;
		
		// Assimilated C/N/P is actually rejected to labile as nematode necromass, to ensure nematode biomass stability
		do reject_to_dom(C_assimilated, N_assimilated, P_assimilated); //dejection rejeté dans le dam
		
		stomack_N <- stomack_N - N_assimilated;
		stomack_P <- stomack_P - P_assimilated;
		
		// Excess N/P is rejected in DIM
		current_pore.dam.dim[0] <- current_pore.dam.dim[0] + stomack_N;
		current_pore.dam.dim[1] <- current_pore.dam.dim[1] + stomack_P;
		
		stomack_N <- 0.0;
		stomack_P <- 0.0;		
	}
	
	
	action reject_to_dom(float C_rejected, float N_rejected, float P_rejected) 
	{
		ask current_pore.dam {
			dom[0] <- dom[0] + N_rejected;
			dom[1] <- dom[1] + P_rejected;
			dom[2] <- dom[2] + C_rejected;
		}
	}
	
 	
 	
 	aspect red_dot {
 		if(awake) {
 			draw cube(0.005) at: location color: #red;
 		} else {
 			draw cube(0.005) at: location color: #brown;
 		}
 	}
}
