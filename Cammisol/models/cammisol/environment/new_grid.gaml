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


	geometry shape <- cube(soil_size);
	float z_bas;
	float taille_min;

	// --- Bilan matière, repris de cammisol_grid.gaml (TODO : vérifier en volume 3D) ---
	float bulk_density <- 1.17#gram / (#cm * #cm * #cm);
	float soil_weight <- bulk_density * ((soil_size * 100) ^ 3) #cm*#cm*#cm;

	// TODO Est ce que ces concentrations sont les memes a toutes les échelles ?

	float C_concentration_in_pom <- 0.02246#gram / (#cm * #cm * #cm);
	float N_concentration_in_pom <- 0.00137#gram / (#cm * #cm * #cm);
	float P_concentration_in_pom <- 0.00070#gram / (#cm * #cm * #cm);

	float C_concentration_in_dom <- (729.0 * 10 ^ -6) #gram / #gram;
	float N_concentration_in_dom <- (60.0 * 10 ^ -6) #gram / #gram;
	float P_concentration_in_dom <- (400.0 * 10 ^ -6) #gram / #gram;

	float N_concentration_in_dim <- (10.1 * 10 ^ -6) #gram / #gram;
	float P_concentration_in_dim <- (10.1 * 10 ^ -6) #gram / #gram;

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

		float C_in_each_dom <- soil_weight * C_concentration_in_dom / pores_count; // TODO Est ce que ces concentrations sont les memes a toutes les échelles ?
		float N_in_each_dom <- soil_weight * N_concentration_in_dom / pores_count; // plutot que en count il faudrait travailler plutot en volume ?
		float P_in_each_dom <- soil_weight * P_concentration_in_dom / pores_count;
		float N_in_each_dim <- soil_weight * N_concentration_in_dim / pores_count;
		float P_in_each_dim <- soil_weight * P_concentration_in_dim / pores_count;

		ask pore_cells {
			PoreParticle pp <- PoreParticle(domain_particle);

			// Voisins réels de la cellule, hors elle-même : le skill APSF renvoie 27 positions

			list<particle> vois <- [];
			loop v over: my_neighbors {
				if (v != self) { add v to: vois; } //on en veut pas que la cellule courante soit voisine
			}

			ask pp {
				add self to: pore_neighbors;   // le pore est son propre premier voisin (comme en 2D)
				loop v over: vois {
					if (v.my_type = PORE) { add PoreParticle(v.domain_particle) to: pore_neighbors; }
				}

				loop v over: vois {
					if (v.my_type = ORGANIC) { add OrganicParticle(v.domain_particle) to: accessible_organics; }
				}

				ask dam {
					dom <- [N_in_each_dom, P_in_each_dom, C_in_each_dom];
					dim <- [N_in_each_dim, P_in_each_dim];
				}
			}
		}

		do init_POM(C: C_concentration_in_pom, N: N_concentration_in_pom, P: P_concentration_in_pom);
	}

	action init_POM(float C, float N, float P) {
		write "Soil C (pom): " + C / (#kg / #m3) + "kg/m3";
		write "Soil N (pom): " + N / (#kg / #m3) + "kg/m3";
		write "Soil P (pom): " + P / (#kg / #m3) + "kg/m3";

		// TODO remplacer par le vrai volume de la cellule APSF (particle.size^3)
		float cell_volume <- (soil_size / division) ^ 3;

		float C_labile_in_each_organic <- labile_rate_pom * C * cell_volume;
		float N_labile_in_each_organic <- labile_rate_pom * N * cell_volume;
		float P_labile_in_each_organic <- labile_rate_pom * P * cell_volume;
		float C_recal_in_each_organic <- (1 - labile_rate_pom) * C * cell_volume;
		float N_recal_in_each_organic <- (1 - labile_rate_pom) * N * cell_volume;
		float P_recal_in_each_organic <- (1 - labile_rate_pom) * P * cell_volume;

		ask OrganicParticle {
			if (!in_pore) {
				C_labile <- C_labile_in_each_organic;
				N_labile <- N_labile_in_each_organic;
				P_labile <- P_labile_in_each_organic;
				C_recalcitrant <- C_recal_in_each_organic;
				N_recalcitrant <- N_recal_in_each_organic;
				P_recalcitrant <- P_recal_in_each_organic;
			}
		}
	}
}

species soil skills: [camisole] {

    init {
        location <- world.location;
        z_bas <- location.z - soil_size / 2;

        // On définit les canvas dans l'ordre des dépendances (C1 avant C0)
        list<string> to_define <- copy(canvases.keys);
        list<string> defined <- [];
        int guard <- 0;

        loop while: (!empty(to_define) and guard < 100) {
            guard <- guard + 1;

            loop canvas_name over: copy(to_define) {
                // noms de canvas utilisés comme cellules dans ce canvas
                list<string> deps <- [];
                loop plan over: list(canvases[canvas_name]) {
                    loop row over: list(plan) {
                        loop c over: list(row) {
                            if (string(c) in canvases.keys) { deps <+ string(c); }
                        }
                    }
                }

                if (deps all_match (each in defined)) {
                    do define_canvas(
                        name: canvas_name,
                        cells: list<list<list<string>>>(canvases[canvas_name])
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
            from_canvas: root_canvas,
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
