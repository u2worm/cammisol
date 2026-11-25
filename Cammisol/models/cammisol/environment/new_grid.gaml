/**
 * Name: new_cammisol_grid
 * Remplace la grille GAMA 2D de cammisol_grid.gaml par une grille APSF 3D fractale.
 * Les espèces domaine (MineralParticle, OrganicParticle, PoreParticle, Dam, MicrobePopulation,
 * ...) ne sont pas modifiées. La couche de construction de la grille et des listes de
 * voisinage change.
 *
 * Phase de ce fichier : câblage grille + biologie, à une seule échelle de déploiement.
 * Les cellules décomposables (agglomérats non subdivisés) sont volontairement ignorées.
 * Les formules de bilan matière (soil_weight, cell_volume) sont reprises telles quelles de
 * la version 2D et PAS encore recalibrées en volume 3D :
 */

model new_cammisol_grid

import "pore_particle.gaml"
import "mineral_particle.gaml"

global {
	string ORGANIC <- "organic";
	string MINERAL <- "mineral";
	string PORE <- "pore";

	// --- Grille APSF ----------------------------------------------------------------------------------
	float soil_size <- 0.2;    // côté du cube de sol, en mètres
	int division <- 3;         // cellules par côté d'un canevas (doit correspondre aux canevas définis)
	int niveaux <- 2;          // profondeur de déploiement : une seule échelle pour cette phase //TODO gérer ce qu'il se passe pour les decomposables
	string soil_file <- "../../../canvases/grid4.json"; // fichier json comprenant le caneva fixe
	map<string, list> canvases <- []; // tableau de cellules
	string root_canvas <- "C0"; // racine du caneva

	// [Canevas] Canevas « déroulés » par niveau, réellement envoyés au plugin (voir unroll_canvases).
	// Chaque canevas du fichier est recopié une fois par niveau où il est utilisé, sous le nom
	// <nom>@<niveau> (ex. "C1@2"). Cela permet de traiter différemment le dernier niveau, et de gérer
	// les canevas qui se référencent eux-mêmes (vrai fractal) sans dépendance circulaire.
	map<string, list> canvases_deployed <- [];
	string root_canvas_deployed <- "";      // nom de la racine : root_canvas + "@1"
	string LEVEL_SEPARATOR <- "@";          // séparateur nom/niveau, interdit dans les noms du fichier
	int min_canvas_size <- 3;               // taille minimale d'un canevas (voir validate_canvases)

	// [Voisinage] true : voisinage par faces (6 directions, multi-échelle, symétrique, sans doublon),
	// fourni par la variable my_face_neighbors du plugin (nécessite le plugin recompilé).
	// false : ancien voisinage de Moore 3D (26 voisins, my_neighbors), conservé pour comparaison.
	bool voisinage_faces <- true;


	geometry shape <- cube(soil_size);
	float z_bas;
	float taille_min;

	// --- Bilan matière, repris de cammisol_grid.gaml (TODO : vérifier en volume 3D) ---
	float bulk_density <- 1.17#gram / (#cm * #cm * #cm);
	float soil_weight <- bulk_density * ((soil_size * 100) ^ 3) #cm*#cm*#cm;

	// TODO Est ce que ces concentrations sont les memes a toutes les échelles ?

	// [Initialisation] POM par masse de sol sec (g/g), ODD Tab. 1 (sol noir). Total = conc. x soil_weight.
	float C_concentration_in_pom <- 0.02246 #gram / #gram;
	float N_concentration_in_pom <- 0.00137 #gram / #gram;
	float P_concentration_in_pom <- 0.00070 #gram / #gram;

	float C_concentration_in_dom <- (72.0 * 10 ^ -6) #gram / #gram;
	float N_concentration_in_dom <- (6.00 * 10 ^ -6) #gram / #gram;
	float P_concentration_in_dom <- (40.0 * 10 ^ -6) #gram / #gram;

	float N_concentration_in_dim <- (1.01 * 10 ^ -6) #gram / #gram;
	float P_concentration_in_dim <- (0.119 * 10 ^ -6) #gram / #gram;

	float labile_rate_pom <- 0.2;

	// --- Affichage ----------------------------------------------------------------------------------
	float coupe <- 1.0;
	float opacite <- 1.0;
	bool voir_mineral <- true;
	bool voir_organic <- true;
	bool voir_agglo <- true;
	bool voir_pores <- false;
	string style_pores <- "wireframe" among: ["wireframe", "transparent"];
	bool voir_cadre <- true;
	bool show_nematodes <- true;

	/*init {
		do load_soil_file();	// charge le Json
		create soil number: 1;   // définit les canevas, initialize_soil_from_canvas, materialize
		//do init_grid();          // câble le domaine CAMMISOL sur les agents 'particle' créés ci-dessus
	}*/

	// ============================================================ // Lecture du fichier // ============================================================
	
	action load_soil_file() {
	    json_file f <- json_file(soil_file);
	    map<string, unknown> data <- map<string, unknown>(f.contents);
	
	    division <- int(data["size"]);
	
	    map<string, unknown> raw_canvases <- map<string, unknown>(data["canvases"]);
	    canvases <- [];
	
	    loop canvas_name over: raw_canvases.keys {
	        canvases[canvas_name] <- list(raw_canvases[canvas_name]);
	    }
	
	    write "Fichier : " + soil_file;
	    write "Taille grille : " + division;
	    write "Canvas disponibles : " + canvases.keys;

	    // [Canevas] Vérifications puis déroulage par niveau, avant tout envoi au plugin
	    do validate_canvases();
	    do unroll_canvases();
	}

	// ============================================================ // Vérification des canevas // ============================================================

	/**
	 * [Canevas] Vérifie la cohérence des canevas lus dans le fichier JSON. Toute anomalie arrête la
	 * simulation avec un message explicite (plutôt qu'une erreur obscure plus loin dans le plugin).
	 */
	action validate_canvases() {
		list<string> base_types <- [PORE, ORGANIC, MINERAL];

		// au moins un canevas, et le canevas racine (root_canvas) existe ;
		if (empty(canvases.keys)) { error "[Canevas] Aucun canevas dans le fichier " + soil_file; }
		if (!(root_canvas in canvases.keys)) {
			error "[Canevas] Le canevas racine '" + root_canvas + "' n'existe pas dans " + soil_file
				+ " (canevas disponibles : " + canvases.keys + ")";
		}
		// n >= 3 : avec le tore, un canevas de 2 de côté donnerait le même voisin
		if (division < min_canvas_size) {
			error "[Canevas] Taille de canevas " + division + " < " + min_canvas_size
				+ " : avec le tore, les voisins x-1 et x+1 seraient la même cellule.";
		}

		loop canvas_name over: canvases.keys {
			// aucun nom de canevas ne contient le séparateur "@" ni ne vaut "pore", "organic" ou "mineral" ;
			if (canvas_name contains LEVEL_SEPARATOR) {
				error "[Canevas] Le nom '" + canvas_name + "' contient '" + LEVEL_SEPARATOR + "', réservé au déroulage par niveau.";
			}
			if (canvas_name in base_types) {
				error "[Canevas] Le nom de canevas '" + canvas_name + "' est réservé (type de cellule).";
			}
			// chaque canevas est un cube n x n x n ;
			// n est le même pour tous les canevas et égal au champ "size" du fichier
			list planes <- list(canvases[canvas_name]);
			if (length(planes) != division) {
				error "[Canevas] '" + canvas_name + "' a " + length(planes) + " plans au lieu de " + division
					+ " (tous les canevas doivent avoir la taille 'size' du fichier).";
			}
			loop i from: 0 to: length(planes) - 1 {
				list rows <- list(planes[i]);
				if (length(rows) != division) {
					error "[Canevas] '" + canvas_name + "', plan " + i + " : " + length(rows) + " lignes au lieu de " + division;
				}
				loop j from: 0 to: length(rows) - 1 {
					list cells <- list(rows[j]);
					if (length(cells) != division) {
						error "[Canevas] '" + canvas_name + "', plan " + i + ", ligne " + j + " : "
							+ length(cells) + " cellules au lieu de " + division;
					}
					// chaque cellule vaut "pore", "organic", "mineral" ou le nom d'un canevas
					loop k from: 0 to: length(cells) - 1 {
						string t <- string(cells[k]);
						if (!(t in base_types) and !(t in canvases.keys)) {
							error "[Canevas] '" + canvas_name + "' [" + i + "][" + j + "][" + k + "] : type de cellule inconnu '" + t + "'";
						}
					}
				}
			}
		}
		write "[Canevas] Vérification OK : " + length(canvases.keys) + " canevas de " + division + "x" + division + "x" + division;
	}

	/**
	 * [Canevas] Déroule l'arbre des canevas niveau par niveau, de 1 (racine) à 'niveaux'.
	 *
	 * un même canevas peut être utilisé à plusieurs niveaux (et même se contenir lui-même).
	 * Or les cellules décomposables ne doivent être remplacées qu'au DERNIER niveau. On crée donc une
	 * copie par niveau, nommée <nom>@<niveau> :
	 *  - à un niveau L, une cellule qui référence le canevas Y devient une référence à "Y@(L+1)" ;
	 *
	 * Seuls les canevas atteignables depuis la racine sont déroulés.
	 */
	action unroll_canvases() {
		if (niveaux < 1) { error "[Canevas] 'niveaux' doit être >= 1 (valeur : " + niveaux + ")"; }

		canvases_deployed <- [];
		root_canvas_deployed <- root_canvas + LEVEL_SEPARATOR + 1;

		// Nombre d'exemplaires de chaque canevas source au niveau courant (la racine existe une fois)
		map<string, int> instances_current <- [root_canvas::1];
		int replaced_in_definitions <- 0;   // cases modifiées dans les copies
		int replaced_in_soil <- 0;          // cellules réellement concernées dans le sol
		list<string> summary <- [];         // ex. "C1@2 x15"
		list<string> used_sources <- [];    // canevas du fichier effectivement utilisés

		loop level from: 1 to: niveaux {
			map<string, int> instances_next <- [];
			loop src over: instances_current.keys {
				int nb_instances <- instances_current[src];
				int replaced_here <- 0;
				list<list<list<string>>> new_planes <- [];
				loop plane over: list(canvases[src]) {
					list<list<string>> new_rows <- [];
					loop row over: list(plane) {
						list<string> new_cells <- [];
						loop c over: list(row) {
							string t <- string(c);
							if (t in canvases.keys) {
								// Cellule décomposable
								if (level < niveaux) {
									add t + LEVEL_SEPARATOR + (level + 1) to: new_cells;
									instances_next[t] <- (instances_next contains_key t ? instances_next[t] : 0) + nb_instances;
								} else {
									add PORE to: new_cells;   // dernier niveau : décomposable remplacée par un pore
									replaced_here <- replaced_here + 1;
								}
							} else {
								add t to: new_cells;          // pore, organique ou minéral : inchangé
							}
						}
						add new_cells to: new_rows;
					}
					add new_rows to: new_planes;
				}
				string unrolled_name <- src + LEVEL_SEPARATOR + level;
				canvases_deployed[unrolled_name] <- new_planes;
				summary <+ unrolled_name + " x" + nb_instances;
				if (!(used_sources contains src)) { used_sources <+ src; }
				replaced_in_definitions <- replaced_in_definitions + replaced_here;
				replaced_in_soil <- replaced_in_soil + replaced_here * nb_instances;
			}
			instances_current <- instances_next;
			if (empty(instances_current)) { break; }   // plus aucune décomposable : arbre terminé
		}

		//write "[Canevas] Déroulage sur " + niveaux + " niveau(x) : " + summary;
		// Canevas du fichier jamais atteints depuis la racine (simple information)
		list<string> unused <- canvases.keys where !(used_sources contains each);
		if (!empty(unused)) { write "[Canevas] Canevas non utilisés depuis la racine '" + root_canvas + "' : " + unused; }
		if (replaced_in_soil > 0) {
			string msg <- "[Canevas] " + replaced_in_soil + " cellule(s) décomposable(s) au dernier niveau (" + niveaux
				+ ") remplacée(s) par des pores dans le sol (" + replaced_in_definitions;
			warn msg;
			write msg;
		}
	}
	
	


	/**
	 * nouveau init_grid crée les agents domaine (via
	 * particle.init_particle) puis construit les listes de voisins
	 */
	action init_grid() {
		
		if (length(canvases.keys) = 0) {
	        error "ERREUR : aucun canvas n'a été chargé.";
	        return;
	    }
		
		ask particle {
			do init_particle();
		}

		list<particle> pore_cells <- particle where (each.my_type = PORE);
		int pores_count <- length(pore_cells);

		// [Initialisation] chaque pore reçoit une part du contenu initial proportionnelle à son
		// volume (part_volume), et non plus 1 / nombre de pores. Les totaux du domaine ne changent pas.
		// TODO Est ce que ces concentrations sont les memes a toutes les échelles ?
		float pores_total_volume <- sum(pore_cells collect (each.size ^ 3));
		ask pore_cells {
			ask PoreParticle(domain_particle) { part_volume <- (size ^ 3) / pores_total_volume; }
		}

		// [Initialisation] Totaux du domaine
		float C_dom_total <- soil_weight * C_concentration_in_dom;
		float N_dom_total <- soil_weight * N_concentration_in_dom;
		float P_dom_total <- soil_weight * P_concentration_in_dom;
		float N_dim_total <- soil_weight * N_concentration_in_dim;
		float P_dim_total <- soil_weight * P_concentration_in_dim;

		ask pore_cells {
			PoreParticle pp <- PoreParticle(domain_particle);

			// [Voisinage] Construction de la liste des voisins de la cellule
			list<particle> vois <- [];
			if (voisinage_faces) {
				// Voisinage par faces : le plugin renvoie déjà une liste sans la cellule elle-même et sans
				// doublon. Entre deux échelles, une grosse cellule voit tous les petits voisins collés à ses
				// faces, et chaque petit voit la grosse : la relation est symétrique.
				// (Les surfaces de contact sont disponibles dans my_face_contact_areas, même ordre,
				// mais ne sont pas encore utilisées : ce sera un changement séparé.)
				loop v over: my_face_neighbors {
					add particle(v) to: vois;
				}
			} else {
				// Ancien voisinage de Moore 3D : le skill APSF renvoie 27 positions, cellule comprise
				loop v over: my_neighbors {
					if (v != self) { add particle(v) to: vois; } //on en veut pas que la cellule courante soit voisine
				}
			}

			ask pp {
				add self to: pore_neighbors;   // le pore est son propre premier voisin (comme en 2D)
				loop v over: vois {
					if (v.my_type = PORE) { add PoreParticle(v.domain_particle) to: pore_neighbors; }
				}

				loop v over: vois {
					if (v.my_type = ORGANIC) { add OrganicParticle(v.domain_particle) to: accessible_organics; }
				}

				float w <- part_volume;   // [Initialisation] part du volume des pores
				ask dam {
					dom <- [N_dom_total * w, P_dom_total * w, C_dom_total * w];
					dim <- [N_dim_total * w, P_dim_total * w];
				}
			}
		}

		// [Voisinage] Diagnostic de la structure.
		list<PoreParticle> all_pores <- list<PoreParticle>(pore_cells collect each.domain_particle);
		list<OrganicParticle> reachable <- remove_duplicates(all_pores accumulate each.accessible_organics);
		int nb_isolated_pores <- length(all_pores where (length(each.pore_neighbors) <= 1)); // seul voisin = lui-même
		int nb_pores_without_organic <- length(all_pores where empty(each.accessible_organics));
		int nb_unreachable_organics <- length(OrganicParticle) - length(reachable);
		/*write "[Voisinage] " + (voisinage_faces ? "par faces (6, multi-échelle)" : "Moore 3D (26)") + " :";
		write "  pores voisins par pore (moyenne, hors lui-même) : "
			+ (empty(all_pores) ? 0.0 : mean(all_pores collect (length(each.pore_neighbors) - 1)));
		write "  organiques accessibles par pore (moyenne) : "
			+ (empty(all_pores) ? 0.0 : mean(all_pores collect length(each.accessible_organics)));
		write "  pores isolés : " + nb_isolated_pores + " / " + length(all_pores);
		write "  pores sans organique voisin : " + nb_pores_without_organic + " / " + length(all_pores);
		write "  organiques inaccessibles (aucun pore voisin) : " + nb_unreachable_organics + " / " + length(OrganicParticle);*/

		/* 
		list<float> pore_sizes <- remove_duplicates(all_pores collect each.size);
		write "[Initialisation] Pores : " + length(all_pores) + " pores, " + length(pore_sizes) + " taille(s), "
			+ "volume de pores = " + (100 * pores_total_volume / (soil_size ^ 3)) + " % du domaine, "
			+ "somme des part_volume = " + sum(all_pores collect each.part_volume);
		do diag_total(label: "C DOM", s: sum(all_pores collect each.dam.dom[2]), target: C_dom_total);
		do diag_total(label: "N DOM", s: sum(all_pores collect each.dam.dom[0]), target: N_dom_total);
		do diag_total(label: "P DOM", s: sum(all_pores collect each.dam.dom[1]), target: P_dom_total);
		do diag_total(label: "N DIM", s: sum(all_pores collect each.dam.dim[0]), target: N_dim_total);
		do diag_total(label: "P DIM", s: sum(all_pores collect each.dam.dim[1]), target: P_dim_total);
		do diag_density(label: "C DOM", values: all_pores collect (each.dam.dom[2] / (each.size ^ 3)));
		do diag_density(label: "N DIM", values: all_pores collect (each.dam.dim[0] / (each.size ^ 3)));*/

		do init_POM(C: C_concentration_in_pom, N: N_concentration_in_pom, P: P_concentration_in_pom);
	}


	action diag_total(string label, float s, float target) {
		float rel <- target = 0.0 ? abs(s) : abs(s - target) / target;
		write "  " + label + " = " + (s / #gram) + " g (cible " + (target / #gram) + " g, écart relatif " + rel + ")";
	}

	// [Initialisation] Rapport max/min d'une densité (quantité / volume) entre cellules : 1 = parfaitement uniforme.
	action diag_density(string label, list<float> values) {
		float ratio <- (empty(values) or min(values) <= 0.0) ? -1.0 : max(values) / min(values);
		write "  densité " + label + " max / min entre pores = " + ratio;
	}

	/**
	 * [Initialisation] répartition de la POM initiale au prorata du VOLUME des cellules organiques.
	 * C, N, P : concentrations par masse de sol sec (g/g, Tableau 1 de l'ODD).
	 *  - total du domaine : X_tot = X x soil_weight (indépendant du canevas et du nombre de niveaux) ;
	 *  - chaque organique reçoit X_tot x size^3 / somme des size^3 : la densité (g/m3) est donc la même
	 *    dans toutes les cellules organiques, quelle que soit leur taille.  //TODO a verifier
	 */
	action init_POM(float C, float N, float P) {
		write "Soil C (pom): " + C / (#gram / #kg) + " g/kg de sol";
		write "Soil N (pom): " + N / (#gram / #kg) + " g/kg de sol";
		write "Soil P (pom): " + P / (#gram / #kg) + " g/kg de sol";

		list<OrganicParticle> organics <- OrganicParticle where (!each.in_pore);
		float total_volume <- sum(organics collect (each.size ^ 3));
		if (empty(organics) or total_volume <= 0.0) {
			warn "[Initialisation] Aucune particule organique : la POM initiale n'est pas attribuée.";
			return;
		}

		// [Initialisation] Totaux du domaine (cibles du test)
		float C_total <- C * soil_weight;
		float N_total <- N * soil_weight;
		float P_total <- P * soil_weight;

		ask organics {
			float w <- (size ^ 3) / total_volume;   // part du volume organique total portée par cette cellule
			C_labile <- labile_rate_pom * C_total * w;
			N_labile <- labile_rate_pom * N_total * w;
			P_labile <- labile_rate_pom * P_total * w;
			C_recalcitrant <- (1 - labile_rate_pom) * C_total * w;
			N_recalcitrant <- (1 - labile_rate_pom) * N_total * w;
			P_recalcitrant <- (1 - labile_rate_pom) * P_total * w;
		}

		// [Initialisation]
		//  somme des organiques = cible (écart relatif attendu < 1e-12) ;
		//  densité C / volume identique dans toutes les cellules (rapport max/min attendu = 1).
		float C_sum <- sum(organics collect (each.C_labile + each.C_recalcitrant));
		float N_sum <- sum(organics collect (each.N_labile + each.N_recalcitrant));
		float P_sum <- sum(organics collect (each.P_labile + each.P_recalcitrant));
		list<float> densities <- organics collect ((each.C_labile + each.C_recalcitrant) / (each.size ^ 3));
		list<float> sizes <- remove_duplicates(organics collect each.size);
		/*write "[Initialisation] POM : " + length(organics) + " organiques, " + length(sizes) + " taille(s), "
			+ "volume organique = " + (100 * total_volume / (soil_size ^ 3)) + " % du domaine";
		write "  masse de sol = " + (soil_weight / #kg) + " kg";
		write "  C POM = " + (C_sum / #gram) + " g (cible " + (C_total / #gram) + " g, écart relatif " + (abs(C_sum - C_total) / C_total) + ")";
		write "  N POM = " + (N_sum / #gram) + " g (cible " + (N_total / #gram) + " g, écart relatif " + (abs(N_sum - N_total) / N_total) + ")";
		write "  P POM = " + (P_sum / #gram) + " g (cible " + (P_total / #gram) + " g, écart relatif " + (abs(P_sum - P_total) / P_total) + ")";
		write "  densité C max / min entre organiques = " + (max(densities) / min(densities));*/
	}
}

species soil skills: [camisole] {

    init {
        location <- world.location;
        z_bas <- location.z - soil_size / 2;

        // On définit les canvas dans l'ordre des dépendances (C1 avant C0)
        // [Canevas] On envoie au plugin les canevas déroulés par niveau (canvases_deployed, construits par
        // unroll_canvases dans load_soil_file), et non plus les canevas bruts du fichier.
        list<string> to_define <- copy(canvases_deployed.keys);
        list<string> defined <- [];
        int guard <- 0;

        loop while: (!empty(to_define) and guard < 100) {
            guard <- guard + 1;

            loop canvas_name over: copy(to_define) {
                // noms de canvas utilisés comme cellules dans ce canvas
                list<string> deps <- [];
                loop plan over: list(canvases_deployed[canvas_name]) {
                    loop row over: list(plan) {
                        loop c over: list(row) {
                            if (string(c) in canvases_deployed.keys) { deps <+ string(c); }
                        }
                    }
                }

                if (deps all_match (each in defined)) {
                    do define_canvas(
                        name: canvas_name,
                        cells: list<list<list<string>>>(canvases_deployed[canvas_name])
                    );
                    defined <+ canvas_name;
                    to_define >- canvas_name;
                    write "Canvas défini : " + canvas_name;
                }
            }
        }

        if (!empty(to_define)) {
            error "Dépendance circulaire ou canvas manquant : " + to_define;
        }

        do initialize_soil_from_canvas(
            from_canvas: root_canvas_deployed,   // [Canevas] racine déroulée (ex. "C0@1")
            default_species: particle,
            size: soil_size
        );
        do materialize(depth: niveaux);

        taille_min <- min(particle collect each.size);

        ask world { do init_grid(); }

        write "Cellules créées : " + length(particle);
        loop t over: remove_duplicates(particle collect each.my_type) {
            write "  " + t + " : " + length(particle where (each.my_type = t));
        }
    }
    aspect cadre { 
    	if (voir_cadre) { 
    		draw cube(soil_size) at: location + {0, 0, -soil_size / 2} wireframe: true color: #red;
    	}
    }
}

species particle skills: [apsf_particle] {
	// Lien vers l'agent domaine CAMMISOL (MineralParticle / OrganicParticle / PoreParticle).
	// Reste nil pour les cellules décomposables (my_type = nom d'un canevas), ignorées pour cette phase.
	agent domain_particle;

	bool est_feuille -> my_type in [MINERAL, ORGANIC, PORE];
	bool est_dessinee -> est_feuille or abs(size - taille_min) < 1e-9;
	point pos_dessin -> location + {0, 0, -size / 2};

	action init_particle() {
		switch my_type {
			match MINERAL {
				create MineralParticle with: (grid_x: i, grid_y: j, grid_z : k,size : size, location: location) { // TODO verifier que k correspond bien a z
					myself.domain_particle <- self;
				}
			}
			match ORGANIC {
				create OrganicParticle with: (grid_x: i, grid_y: j,grid_z : k, size : size, location: location) {
					myself.domain_particle <- self;
				}
			}
			match PORE {

				Dam pore_dam;


				create Dam { pore_dam <- self; }
				create PoreParticle with: (
					 dam: pore_dam,
					grid_x: i, grid_y: j, grid_z : k, size : size, location: location
				) {
					myself.domain_particle <- self;
				}
			}
			// cellule décomposable (nom de canevas) : pas de domain_particle pour cette phase TODO ajouter les decomposables
		}
	}

	aspect base {
		if (est_dessinee and location.z <= z_bas + coupe * soil_size) {
			int alpha <- int(opacite * 255);
			switch my_type {
				match MINERAL {
					if (voir_mineral) { draw cube(size) at: pos_dessin color: rgb(0, 0, 0, alpha) border: #white; }
				}
				match ORGANIC {
					if (voir_organic) { draw cube(size) at: pos_dessin color: rgb(0, 160, 60, alpha) border: #white; }
				}
				match PORE {
					if (voir_pores) {
						if (style_pores = "wireframe") {
							draw cube(size) at: pos_dessin wireframe: true color: #darkblue;
						} else {
							draw cube(size) at: pos_dessin color: rgb(80, 120, 255, int(opacite * 60)) border: #darkblue;
						}
					}
				}
				default {
					if (voir_agglo) { draw cube(size) at: pos_dessin color: rgb(128, 128, 128, alpha) border: #white; }
				}
			}
		}
	}
}

experiment grille type: gui {
	parameter "Coupe (fraction de hauteur)" var: coupe min: 0.05 max: 1.0 step: 0.05 category: "Affichage";
	parameter "Opacité" var: opacite min: 0.05 max: 1.0 step: 0.05 category: "Affichage";
	parameter "Afficher minéral" var: voir_mineral category: "Affichage";
	parameter "Afficher organique" var: voir_organic category: "Affichage";
	parameter "Afficher décomposable" var: voir_agglo category: "Affichage";
	parameter "Afficher pores" var: voir_pores category: "Affichage";
	parameter "Style des pores" var: style_pores category: "Affichage";
	parameter "Afficher le cadre" var: voir_cadre category: "Affichage";

	output {
		display vue type: opengl background: #white {
			species soil aspect: cadre;
			species particle aspect: base;
		}
	}
}
