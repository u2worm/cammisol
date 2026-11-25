/**
 * Name: nematode
 * Nématodes bactérivores : ils se déplacent entre pores, prédatent les populations microbiennes,
 * respirent et minéralisent le N et le P en excès.
 * La biomasse des nématodes est constante (ni reproduction ni mort) : ce qu'un nématode assimile
 * est aussitôt rejeté sous forme de matière organique labile (état stationnaire).
 * Author: pbreugno
 */

model nematode

import "../microbes/microbes.gaml"

global {
	float nematode_CO2_emissions <- 0.0#gram;                     // CO2 cumulé émis par tous les nématodes
	float nematode_predation_rate <- 3.0e-4 #gram/#day;      // C prédaté par nématode : ~100k bactéries (1e-12 g) par jour
	float nematode_C_N <- 6.0;                                    // rapport C:N de la biomasse des nématodes
	float nematode_C_P <- 47.0;                                   // rapport C:P de la biomasse des nématodes
	float nematode_CUE <- 0.16;                                   // part du C de l'estomac conservée ; le reste est respiré
	float nematode_sleep_threshold <- 0.1;                        // dort si le C perçu < 10 % de la demande (~10k bactéries/j, réf. Colman)
	float nematode_rejected_microbe_necromass_rate <- 0.5;        // part du prédaté rejetée directement en DOM (nécromasse microbienne)
	// Taille (côté) minimale d'un pore accessible aux nématodes, pour le placement initial et les
	// déplacements. 0 = aucune exclusion.
	float taille_min_pore_nematode <- 0.0 #m;
}

species Nematode schedules:[]
{
	// Contenu de l'estomac (matière prédatée en cours de digestion)
	float stomack_C;
	float stomack_N;
	float stomack_P;
	
	float requested_C;     // C demandé pour le pas de temps
	float perceived_C;     // C microbien (structural + cytosol) présent dans le pore courant
	
	PoreParticle current_pore;
	
	bool awake <- true;
	
	/**
	 * Cycle d'un pas : déplacement, perception, puis, si le nématode est éveillé,
	 * prédation, respiration et anabolisme.
	 */
	action life() {
	    requested_C <- nematode_predation_rate * local_step;
	    
	    do move();
	    perceived_C <- perceive(); // après le déplacement
	    
	    awake <- requested_C > 0 and (perceived_C / requested_C) >= nematode_sleep_threshold;
	    if awake { do predation(); do respirate(); do anabolize(); }
	}
	
	/** C microbien total (structural + cytosol, populations actives et dormantes) du pore courant. */
	float perceive(){
		return sum(current_pore.populations collect(each.C + each.cytosol_C));
	}
	
	/**
	 * Se déplace vers le pore voisin (ou reste dans le pore courant) qui contient le plus de C
	 * microbien. Les pores plus petits que taille_min_pore_nematode sont ignorés.
	 */
	action move()
	{
		list<PoreParticle> neighbors_pores <- current_pore.pore_neighbors where (each.size >= taille_min_pore_nematode);
		
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
	}
	
	/**
	 * Prédation : prélève min(C perçu, C demandé) dans les populations du pore.
	 * Prédation proportionnelle : chaque population (et, en son sein, le structural et le cytosol)
	 * perd la même fraction total_C_predated / perceived_C de son C, avec le N et le P correspondants.
	 * Une part nematode_rejected_microbe_necromass_rate est rejetée directement en DOM ; le reste
	 * va dans l'estomac.
	 */
	action predation()
 	{
		float total_C_predated <- min([perceived_C, requested_C]);
		float total_N_predated <- 0.0;
		float total_P_predated <- 0.0;
		
		ask current_pore.populations
		{	
			float total_C_predated_from_this_population <- (total_C_predated / myself.perceived_C) * (C + cytosol_C);
			
			if (total_C_predated_from_this_population > 0.0) {
				// Part prise sur la biomasse structurale (N et P au rapport de la biomasse)
				float C_predated_from_structural <- C/(C+cytosol_C) * total_C_predated_from_this_population;
				float N_predated_from_structural <- C > 0.0 ? C_predated_from_structural * N/C : 0.0;   // sécurité division par 0
				float P_predated_from_structural <- C > 0.0 ? C_predated_from_structural * P/C : 0.0;
							
				C <- C - C_predated_from_structural;
				N <- N - N_predated_from_structural;
				P <- P - P_predated_from_structural;
				
				total_N_predated <- total_N_predated + N_predated_from_structural;
				total_P_predated <- total_P_predated + P_predated_from_structural;
				
				// Part prise sur le cytosol (N et P au rapport du cytosol)
				float C_predated_from_cytosol <- total_C_predated_from_this_population - C_predated_from_structural;
				
				if(C_predated_from_cytosol > 0.0) {
					float N_predated_from_cytosol <- C_predated_from_cytosol * cytosol_N / cytosol_C;	
					float P_predated_from_cytosol <- C_predated_from_cytosol * cytosol_P / cytosol_C;
						
					cytosol_C <- cytosol_C - C_predated_from_cytosol;
					cytosol_N <- cytosol_N - N_predated_from_cytosol;
					cytosol_P <- cytosol_P - P_predated_from_cytosol;
					
					total_N_predated <- total_N_predated + N_predated_from_cytosol;
					total_P_predated <- total_P_predated + P_predated_from_cytosol;
				}
			}
		}
		
		// Nécromasse microbienne (prédation « gaspilleuse ») : rejetée directement en DOM
		float microbe_necromass_C <- total_C_predated*nematode_rejected_microbe_necromass_rate;
		float microbe_necromass_N <- total_N_predated*nematode_rejected_microbe_necromass_rate;
		float microbe_necromass_P <- total_P_predated*nematode_rejected_microbe_necromass_rate;
		
		do reject_to_dom(microbe_necromass_C, microbe_necromass_N, microbe_necromass_P);
		
		stomack_C <- stomack_C + total_C_predated - microbe_necromass_C;
		stomack_N <- stomack_N + total_N_predated - microbe_necromass_N;
		stomack_P <- stomack_P + total_P_predated - microbe_necromass_P;
 	}
	
	/** Respiration : une part (1 - CUE) du C de l'estomac est émise en CO2. */
	action respirate()
	{
		float C_respiration <- (1-nematode_CUE) * stomack_C;
		stomack_C <- stomack_C  - C_respiration;
		nematode_CO2_emissions <- nematode_CO2_emissions + C_respiration;
	}
	
	/**
	 * Anabolisme : le C, N et P de l'estomac utilisables aux rapports C:N et C:P des nématodes
	 * forment de la nouvelle biomasse. La biomasse des nématodes étant constante, cette biomasse est
	 * aussitôt rejetée en MO labile des organiques voisines (reject_to_labile).
	 * Le N et le P en excès sont minéralisés (déjections) dans le DIM. Le C non utilisé reste dans
	 * l'estomac.
	 */
	action anabolize()
	{
		float C_assimilated <- min(stomack_C, stomack_N * nematode_C_N, stomack_P * nematode_C_P);
		float N_assimilated <- C_assimilated/nematode_C_N;
		float P_assimilated <- C_assimilated/nematode_C_P;
		
		stomack_C <- stomack_C - C_assimilated;
		stomack_N <- stomack_N - N_assimilated;
		stomack_P <- stomack_P - P_assimilated;
		
		do reject_to_labile(C_assimilated, N_assimilated, P_assimilated);
		
		// Déjections : N et P en excès minéralisés dans le DIM
		do reject_to_dim(stomack_N, stomack_P);
		
		stomack_N <- 0.0;
		stomack_P <- 0.0;		
	}
	
	/** Ajoute de la matière organique au DOM du pore courant. */
	action reject_to_dom(float C_rejected, float N_rejected, float P_rejected) 
	{
		ask current_pore.dam {
			dom[0] <- dom[0] + N_rejected;
			dom[1] <- dom[1] + P_rejected;
			dom[2] <- dom[2] + C_rejected;
		}
	}
	
	/** Ajoute de la matière minérale (N, P) au DIM du pore courant. */
	action reject_to_dim(float N_rejected, float P_rejected) 
	{
		ask current_pore.dam {
			dim[0] <- dim[0] + N_rejected;
			dim[1] <- dim[1] + P_rejected;
		}
	}
	
	/**
	 * Ajoute de la matière organique à la POM labile des particules organiques accessibles depuis le
	 * pore courant, au prorata de leur volume (size^3). Si le pore n'en touche aucune, la matière
	 * va dans le DOM du pore.
	 */
	action reject_to_labile(float C_rejected, float N_rejected, float P_rejected)
	{
		list<OrganicParticle> organics <- current_pore.accessible_organics;
		float total_volume <- sum(organics collect (each.size ^ 3));
		if (empty(organics) or total_volume <= 0.0) {
			do reject_to_dom(C_rejected, N_rejected, P_rejected);
			return;
		}
		ask organics {
			float w <- (size ^ 3) / total_volume;
			C_labile <- C_labile + C_rejected * w;
			N_labile <- N_labile + N_rejected * w;
			P_labile <- P_labile + P_rejected * w;
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
