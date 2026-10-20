/**
* Name: camisol
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model cammisol

import "environment/new_grid.gaml"
import "nematode/nematode.gaml"

global {
	bool do_opti <- true; // temporaire
	int nematodes_count <- 0;
		
	// 5E8 -> 5E9 bacterie / gramme de sol
	/*
	 * A modifier................................................ 
	 */
	 // TODO: poids total de bactérie par gramme de sol * surface du modèle
	// biomasse : Table 2 de l'ODD (sol noir : 124 mg C/kg)
	float microbles_C_concentration <- 124e-6 #gram / #gram;
	float total_initial_bacteria_weight <- microbles_C_concentration * soil_weight;
	
	// TODO: what rates?
	float init_O_rate <- 0.1;
	float init_F_rate <- 0.2;
	float init_M_rate <- 0.7;
	
	float rain_diffusion_rate <- 0.1;
	float rain_period <- 7#days;
	 
	init {
		do load_soil_file();
		create soil number: 1; // définit les canevas, initialize_soil_from_canvas, materialize pour APSF
		
		//do init_grid(); pas utile
		do init_enzymatic_optimisation();
		do init_enzymes();
		// Counts the number of PORES after the initialization of the grid
		int pores_count <- length(PoreParticle);
		int mineral_count <- length(MineralParticle);
		int organic_count <- length(OrganicParticle);
		

		ask PoreParticle {
			// The carrying capacity of each pore is equal to 10 times the initial bacteria population
			carrying_capacity <- 10 * total_initial_bacteria_weight / pores_count;
			
			create O_Strategist with: [C::init_O_rate * total_initial_bacteria_weight / length(PoreParticle)]{
				add self to:myself.populations;
			}
			create F_Strategist with: [C::init_F_rate * total_initial_bacteria_weight / length(PoreParticle)]{
				add self to:myself.populations;
			}
			create M_Strategist with: [C::init_M_rate * total_initial_bacteria_weight / length(PoreParticle)]{
				add self to:myself.populations;
			}
		}
		
		create Nematode number: nematodes_count
		{
			self.current_pore <- one_of(PoreParticle); 
			location <- current_pore.location;
		}
		if (do_opti){
			ask PoreParticle {
				ask populations {
					do update();
					write "Initial enzymes optimization for " + self;
					do optimize_enzymes(myself.dam, myself.accessible_organics);
				}
			}
		}
		write "number of PoreParticle : " + pores_count;
		write "number of MineralParticle : " + mineral_count;
		write "number of OrganicParticle : " + organic_count;

	}
	


	reflex {
		ask shuffle(Nematode) {
			do life();
		}
		ask shuffle(PoreParticle) {
			ask populations {
				if flip(local_step / enzymes_optimization_period) {
					do update();
					do optimize_enzymes(myself.dam, myself.accessible_organics);
				}
			}
			do decompose();
			do microbe_life();
		}
	}
}

experiment base_cammisol_output {
	parameter "Soil size" category: "Environment" var:soil_size;
	parameter "Nematodes count" category: "Environment" var:nematodes_count;
	//parameter "Organic particle rate" category: "Environment" var:organic_rate;
	//parameter "Mineral particle rate" category: "Environment" var:mineral_rate;
	parameter "Init C (FOM)" category: "Environment" var:C_concentration_in_pom;
	parameter "Init N (FOM)" category: "Environment" var:N_concentration_in_pom;
	parameter "Init P (FOM)" category: "Environment" var:P_concentration_in_pom;
	parameter "Init labile rate (FOM)" category: "Environment" var:labile_rate_pom;
	
	parameter "Nematodes predation rate" category: "Nematode" var:nematode_predation_rate;
	parameter "Nematodes CUE" category: "Nematode" var:nematode_CUE;
	parameter "Nematodes C/N" category: "Nematode" var:nematode_C_N;
	parameter "Nematodes C/P" category: "Nematode" var:nematode_C_P;
	
	parameter "Dividing time (O)" category: "O Strategists" var:dividing_time_O;
	parameter "CUE (O)" category: "O Strategists" var:carbon_use_efficiency_O;
	parameter "Minimum active rate (O)" category: "O Strategists" var:minimum_active_rate_O;
	
	parameter "Dividing time (F)" category: "F Strategists" var:dividing_time_F;
	parameter "CUE (F)" category: "F Strategists" var:carbon_use_efficiency_F;
	parameter "Minimum active rate (F)" category: "F Strategists" var:minimum_active_rate_F;
	
	parameter "Dividing time (M)" category: "M Strategists" var:dividing_time_M;
	parameter "CUE (M)" category: "M Strategists" var:carbon_use_efficiency_M;
	parameter "Minimum active rate (M)" category: "M Strategists" var:minimum_active_rate_M;
	
	// parametres de la vue 3D
	parameter "Coupe (fraction de hauteur)" var: coupe min: 0.05 max: 1.0 step: 0.05 category: "Affichage";
	parameter "Opacité" var: opacite min: 0.05 max: 1.0 step: 0.05 category: "Affichage";
	parameter "Afficher minéral" var: voir_mineral category: "Affichage";
	parameter "Afficher organique" var: voir_organic category: "Affichage";
	parameter "Afficher décomposable" var: voir_agglo category: "Affichage";
	parameter "Afficher pores" var: voir_pores category: "Affichage";
	parameter "Style des pores" var: style_pores category: "Affichage";
	parameter "Afficher le cadre" var: voir_cadre category: "Affichage";
	parameter "Show nematodes" var: show_nematodes category: "Affichage";
	
	reflex update_particle_color {
		float max_population <- 0.0;
		ask PoreParticle {
			float population <- sum(populations collect (each.C + each.cytosol_C));
			if(population > max_population) {
				max_population <- population;
			}
		}
		if(max_population > 0.0) {
			ask particle {
				if(my_type = PORE) {
					color <- rgb(0, 0, 255 * sum(PoreParticle(domain_particle).populations collect (each.C + each.cytosol_C))/max_population);
				}
			}
		}
	}
	
	output {
		/*display vue type: opengl background: #white {
			species soil aspect: cadre;
			species particle aspect: base;
			species Nematode aspect: red_dot visible: show_nematodes;
		}*/
		
		display "CO2" type: java2D {
			chart "Co2" type:series {
				data "microbe_CO2_emissions (g)" value: microbe_CO2_emissions style:spline; 
				data "nematode_CO2_emissions (g)" value: nematode_CO2_emissions style:spline;
			}
		}
		

		display "dom" type: java2D {
			chart "dom composition" type:series {
				data "N dom (g)" value: sum(Dam collect each.dom[0])/#gram style:spline marker:false thickness:3;
				data "P dom (g)" value: sum(Dam collect each.dom[1])/#gram style:spline marker:false thickness:3;
				data "C dom (g)" value: sum(Dam collect each.dom[2])/#gram style:spline marker:false thickness:3;
			}
		}
		
		display "dim" type: java2D {
			chart "dim composition" type:series {
				data "N dim (g)" value: sum(Dam collect each.dim[0])/#gram style:spline marker:false thickness:3;
				data "P dim (g)" value: sum(Dam collect each.dim[1])/#gram style:spline marker:false thickness:3;
			}
		}
		
		display "organics_labile" type:java2D {
			chart "Organics labile" type:series {
				data "C labile (g)" value: sum(OrganicParticle collect each.C_labile)/#gram style:spline marker:false thickness:3;
				data "N labile (g)" value: sum(OrganicParticle collect each.N_labile)/#gram style:spline marker:false thickness:3;
				data "P labile (g)" value: sum(OrganicParticle collect each.P_labile)/#gram style:spline marker:false thickness:3;
			}
		}
		
		display "organics_recalcitrant" type:java2D {
			chart "Organics recalcitrant" type:series {
				data "C recalcitrant (g)" value: sum(OrganicParticle collect each.C_recalcitrant)/#gram style:spline marker:false thickness:3;
				data "N recalcitrant (g)" value: sum(OrganicParticle collect each.N_recalcitrant)/#gram style:spline marker:false thickness:3;
				data "P recalcitrant (g)" value: sum(OrganicParticle collect each.P_recalcitrant)/#gram style:spline marker:false thickness:3;
			}
		}
		
		
		display "populations" type:java2D {
			chart "Bacteria populations" type:series {
				data "O (g)" value: sum(O_Strategist collect each.C)/#gram style:spline color: #red marker:false thickness:3;
				data "F (g)" value: sum(F_Strategist collect each.C)/#gram style:spline color: #green marker:false thickness:3;
				data "M (g)" value: sum(M_Strategist collect each.C)/#gram style:spline color: #blue marker:false thickness:3;
			}
		}
	}
}

experiment test_microbes parent:base_cammisol_output {
	float C_soil <- 1#gram/(#cm#cm);
	float N_soil <- 0.1#gram/(#cm#cm);
	float P_soil <- 0.05#gram/(#cm#cm);
	
	//parameter "Grid size" var:grid_size;
	parameter "Nematodes" var:nematodes_count;
	parameter "C soil" var:C_soil;
	parameter "N soil" var:N_soil;
	parameter "P soil" var:P_soil;
	parameter "Initial O population rate" var:init_O_rate;
	parameter "Initial F population rate" var:init_F_rate;
	parameter "Initial M population rate" var:init_M_rate;
	
	parameter "Enzymes optimization period" var:enzymes_optimization_period init:local_step;
	
	init {
		ask simulation {
			do init_POM(myself.C_soil, myself.N_soil, myself.P_soil);		
		}
	}
}


experiment cammisol parent:base_cammisol_output {
	// Test to reactivate bacterias
	reflex add_N_P when: local_cycle mod 300 = 0 {
		ask simulation {
			// Redistributes the initial quantity of C/N/P in all pores.
			// Notice that this has no scientific meaning, and is used only for test purpose.
			// This might however represent a kind of fertilizer added to the soil.
			int pores_count <- length(PoreParticle);
			float carbon_in_all_pore <- soil_weight * C_concentration_in_dom; 
			float azote_in_all_pore <- soil_weight * N_concentration_in_dom;
			float phosphore_in_all_pore <- soil_weight * P_concentration_in_dom;
			
			ask PoreParticle {
				ask dam {
					float carbon_in_pore <- carbon_in_all_pore / pores_count;
					float azote_in_pore <- azote_in_all_pore / pores_count;
					float phosphore_in_pore <- phosphore_in_all_pore / pores_count;
					
					dom[0] <- dom[0] + azote_in_pore;
					dom[1] <- dom[1] + phosphore_in_pore;
					dom[2] <- dom[2] + carbon_in_pore;
				}
			}
		}
	}
}


//	reflex rain when: rnd(1.0) < (1/rain_period)*local_step {
////		write "It's raining today!";
//		map<species<MicrobePopulation>, MicrobePopulation> extracted_populations;
//		loop bacteria_type over: bacteria_types {
//			create MicrobePopulation with: (
//				C: 0.0, cytosol_C: 0.0,
//				N: 0.0, cytosol_N: 0.0,
//				P: 0.0, cytosol_P: 0.0
//			) {
//				extracted_populations[bacteria_type] <- self;
//			}
//		}
//		
//		ask PoreParticle {
//			ask populations {
//				ask extracted_populations[species(self)] {
//					float extracted_C <- myself.C * rain_diffusion_rate;
//					self.C <- self.C + extracted_C;
//					myself.C <- myself.C - extracted_C;
//					
//					float extracted_P <- myself.P * rain_diffusion_rate;
//					self.P <- self.P + extracted_P;
//					myself.P <- myself.P - extracted_P;
//					
//					float extracted_N <- myself.N * rain_diffusion_rate;
//					self.N <- self.N + extracted_N;
//					myself.N <- myself.N - extracted_N;
//					
//					float extracted_cytosol_C <- myself.cytosol_C * rain_diffusion_rate;
//					self.cytosol_C <- self.cytosol_C + extracted_cytosol_C;
//					myself.cytosol_C <- myself.cytosol_C - extracted_cytosol_C;
//					
//					float extracted_cytosol_P <- myself.cytosol_P * rain_diffusion_rate;
//					self.cytosol_P <- self.cytosol_P + extracted_cytosol_P;
//					myself.cytosol_P <- myself.cytosol_P - extracted_cytosol_P;
//					
//					float extracted_cytosol_N <- myself.cytosol_N * rain_diffusion_rate;
//					self.cytosol_N <- self.cytosol_N + extracted_cytosol_N;
//					myself.cytosol_N <- myself.cytosol_N - extracted_cytosol_N;
//				}
//			}
//		}
//		
//		ask PoreParticle {
//			ask populations {
//				MicrobePopulation extracted_pop <- extracted_populations[species(self)];
//				self.C <- self.C + extracted_pop.C/length(PoreParticle);
//				self.N <- self.N + extracted_pop.N/length(PoreParticle);
//				self.P <- self.P + extracted_pop.P/length(PoreParticle);
//				self.cytosol_C <- self.cytosol_C + extracted_pop.cytosol_C/length(PoreParticle);
//				self.cytosol_N <- self.cytosol_N + extracted_pop.cytosol_N/length(PoreParticle);
//				self.cytosol_P <- self.cytosol_P + extracted_pop.cytosol_P/length(PoreParticle);
//			}
//		}
//	
//		loop bacteria_type over: bacteria_types {
//			ask extracted_populations[bacteria_type] {
//				do die;
//			}
//		}
//	}

