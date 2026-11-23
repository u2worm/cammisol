/**
 * Name: soil_factory_bench
 * Banc d'essai du générateur de sol SoilFactory (plugin camisole), sans biologie.
 * Génère un sol avec initialize_soil, puis vérifie et mesure, niveau par niveau :
 *  - les volumes visés par le plugin, comparés à ceux recalculés ici à partir des mêmes classes ;
 *  - les volumes obtenus par le tirage retenu, comparés aux volumes visés ;
 *  - les cellules créées par materialize, comparées au contenu du premier canevas.
 * Fichiers écrits dans Cammisol/bench_results/soil_factory/ :
 *  - soil_factory_runs.csv : une ligne par run ;
 *  - <run>_niveaux.csv     : une ligne par niveau.
 * Volumes en fraction du volume du cube. Nécessite le plugin recompilé (actions soil_factory_report
 * et soil_factory_rmse).
 */

model soil_factory_bench

global {
	// [Bench] Sol généré
	float soil_size <- 0.2;            // côté du cube (m)
	int division <- 10;                // cellules par côté d'un canevas
	float densite_apparente <- 1.17;   // [Initialisation] g/cm3 (ODD), passée SANS unité GAMA au plugin
	int nb_tirages <- 20000;           // number_of_try : nombre de tirages comparés par SoilFactory
	int profondeur <- 1;               // niveaux créés par materialize
	float tolerance <- 1e-6;           // écart toléré (fraction du cube) pour les tests

	// [Bench] Constantes du plugin (SoilConstant), reprises pour le calcul indépendant
	float densite_minerale <- 2.62;    // g/cm3
	float densite_mo <- 0.8;           // g/cm3
	float facteur_mo <- 1.72;          // MO = 1,72 x C organique

	// [Initialisation] Granulométrie du fichier de test du plugin (TestSoilData), origine à confirmer.
	// Une classe = [masse en g/kg de sol, borne min, borne max] ; borne max 0 = sans borne (jusqu'au cube).
	list<list<float>> classes_minerales <- [
		[164.9, 0.2 #mm, 2 #mm],        // sables grossiers
		[92.5, 0.05 #mm, 0.2 #mm],      // sables fins (bornes supposées : classe en double dans TestSoilData)
		[86.5, 0.02 #mm, 0.05 #mm],     // limons grossiers
		[232.1, 0.002 #mm, 0.02 #mm],   // limons fins
		[423.9, 0.0, 0.002 #mm]         // argiles
	];
	// [Initialisation] Matière organique : masse en g de C organique par kg de sol
	list<list<float>> classes_organiques <- [
		[0.5, 2 #mm, 0.0],              // classe sans borne max
		[0.04, 0.2 #mm, 2 #mm],
		[2.03, 0.05 #mm, 0.2 #mm],
		[28.13, 0.0, 0.05 #mm]
	];

	// [Bench] Sorties
	string bench_dossier <- "../../../bench_results/soil_factory/";
	string bench_etiquette <- "";      // ajoutée au nom des fichiers, vide par défaut
	string bench_graine;
	string bench_prefixe;

	geometry shape <- cube(soil_size);

	init {
		bench_graine <- replace(string(seed), ".", "");
		bench_prefixe <- "sf_t" + nb_tirages + "_n" + profondeur
			+ (bench_etiquette = "" ? "" : "_" + bench_etiquette) + "_s" + bench_graine;
		create generateur number: 1;
	}

	// [Bench] Part (g/kg) des classes comprise entre lo et hi (cm), masse répartie uniformément selon le
	// diamètre, classe sans borne max limitée au cube : même règle que le plugin, écrite indépendamment.
	float part_classes(list<list<float>> classes, float lo, float hi, float cube_cm) {
		float s <- 0.0;
		loop c over: classes {
			float cmin <- c[1] / #cm;
			float cmax <- c[2] = 0.0 ? cube_cm : c[2] / #cm;
			float haut <- min([hi, cmax]);
			float bas <- max([lo, cmin]);
			if (cmax > cmin and haut > bas) { s <- s + c[0] * (haut - bas) / (cmax - cmin); }
		}
		return s;
	}

	// [Bench] Nombre de niveaux attendu : on divise le cube jusqu'à la plus petite borne max des classes
	int niveaux_attendus(float cube_cm) {
		list<float> bornes <- [];
		loop c over: classes_minerales + classes_organiques {
			if (c[2] > 0.0) { bornes <+ c[2] / #cm; }
		}
		float borne <- min(bornes);
		int n <- 0;
		float t <- cube_cm;
		loop while: true {
			t <- t / division;
			n <- n + 1;
			if (t <= borne * (1 + 1e-9)) { break; }
		}
		return n + 1;
	}

	action analyser(generateur g) {
		float cube_cm <- soil_size / #cm;
		float vol_cube <- cube_cm ^ 3;                           // cm3
		float masse_sol <- vol_cube * densite_apparente;         // g
		float f_min <- masse_sol / 1000 / densite_minerale / vol_cube;           // (g/kg) -> fraction du cube
		float f_mo <- masse_sol / 1000 * facteur_mo / densite_mo / vol_cube;     // (g C/kg) -> fraction du cube
		int nb <- length(g.rapport);

		// [Bench] Test 1 : nombre de niveaux
		int nb_attendu <- niveaux_attendus(cube_cm);
		bool test_niveaux <- nb = nb_attendu;

		// [Bench] Test 2 : volumes visés par niveau = volumes recalculés ici
		// [Bench] Test 3 : somme des volumes visés = volume total des classes (aucune masse perdue ni comptée deux fois)
		float ecart_vise <- 0.0;
		float min_vise_total <- 0.0;
		float mo_visee_total <- 0.0;
		list<list<float>> attendus <- [];
		loop r over: g.rapport {
			float lo <- r[1] / #cm;
			float hi <- min([r[2] / #cm, cube_cm]);
			float min_att <- part_classes(classes_minerales, lo, hi, cube_cm) * f_min;
			float mo_att <- part_classes(classes_organiques, lo, hi, cube_cm) * f_mo;
			attendus <+ [min_att, mo_att];
			ecart_vise <- max([ecart_vise, abs(r[3] - min_att), abs(r[4] - mo_att)]);
			min_vise_total <- min_vise_total + r[3];
			mo_visee_total <- mo_visee_total + r[4];
		}
		bool test_vise <- ecart_vise < tolerance;
		float min_total <- sum(classes_minerales collect each[0]) * f_min;
		float mo_total <- sum(classes_organiques collect each[0]) * f_mo;
		float ecart_conservation <- max([abs(min_vise_total - min_total), abs(mo_visee_total - mo_total)]);
		bool test_conservation <- ecart_conservation < tolerance;

		// [Bench] Test 4 : cellules du premier niveau créées par materialize = contenu du canevas racine
		float taille_1 <- soil_size / division;
		// [Bench] le plugin renvoie size en simple précision (float Java) : tolérance relative de 1e-5
		list<cellule> niveau_1 <- cellule where (abs(each.size - taille_1) < 1e-5 * taille_1);
		list<int> compte_mat <- [
			length(niveau_1 where (each.my_type = "mineral")),
			length(niveau_1 where (each.my_type = "organic")),
			length(niveau_1 where (each.my_type = "pore")),
			length(niveau_1 where !(each.my_type in ["mineral", "organic", "pore"]))
		];
		list<int> compte_canevas <- [int(g.rapport[0][8]), int(g.rapport[0][9]), int(g.rapport[0][10]), int(g.rapport[0][11])];
		bool test_materialize <- compte_mat = compte_canevas;

		// [Bench] Mesures : porosité visée (1 - solides visés), pores obtenus, volume non résolu sous le dernier niveau
		float porosite_visee <- 1 - min_vise_total - mo_visee_total;
		float porosite_obtenue <- sum(g.rapport collect each[7]);
		float reste <- 1 - sum(g.rapport collect (each[5] + each[6] + each[7]));

		write "========== SoilFactory : " + bench_prefixe + " ==========";
		write "Génération : " + (g.duree_ms / 1000) + " s pour " + nb_tirages + " tirages";
		write "[Test] niveaux : " + nb + " (attendu " + nb_attendu + ") -> " + (test_niveaux ? "OK" : "ÉCHEC");
		write "[Test] volumes visés = recalculés : écart max " + ecart_vise + " -> " + (test_vise ? "OK" : "ÉCHEC");
		write "[Test] conservation des classes : écart " + ecart_conservation + " -> " + (test_conservation ? "OK" : "ÉCHEC");
		write "[Test] materialize = canevas racine (minéral, MO, pore, décomposable) : " + compte_mat + " / "
			+ compte_canevas + " -> " + (test_materialize ? "OK" : "ÉCHEC");
		write "Niveau | cellule | minéral visé / obtenu | MO visée / obtenue | pores obtenus | canevas (min, MO, pore, déc.)";
		loop r over: g.rapport {
			int i <- int(r[0]);
			write "C" + i + " | " + (r[1] = 0.0 ? "<= " + (r[2] / #mm * 1000) : string(r[1] / #mm * 1000)) + " µm | "
				+ ((100 * r[3]) with_precision 3) + " / " + ((100 * r[5]) with_precision 3) + " % | "
				+ ((100 * r[4]) with_precision 3) + " / " + ((100 * r[6]) with_precision 3) + " % | "
				+ ((100 * r[7]) with_precision 3) + " % | "
				+ [int(r[8]), int(r[9]), int(r[10]), int(r[11])];
		}
		write "Porosité visée " + ((100 * porosite_visee) with_precision 2) + " %, pores obtenus "
			+ ((100 * porosite_obtenue) with_precision 2) + " %, non résolu sous le dernier niveau "
			+ ((100 * reste) with_precision 3) + " %, RMS " + ((100 * g.rmse) with_precision 3) + " % du cube";

		// [Bench] Export
		string f_runs <- bench_dossier + "soil_factory_runs.csv";
		if (!file_exists(f_runs)) {
			save ["run", "graine", "nb_tirages", "profondeur", "niveaux", "duree_s", "rmse", "porosite_visee",
				"pores_obtenus", "non_resolu", "test_niveaux", "test_vise", "test_conservation", "test_materialize"]
				to: f_runs format: "csv" rewrite: true header: false;
		}
		save [bench_prefixe, bench_graine, nb_tirages, profondeur, nb, g.duree_ms / 1000, g.rmse, porosite_visee,
			porosite_obtenue, reste, test_niveaux, test_vise, test_conservation, test_materialize]
			to: f_runs format: "csv" rewrite: false header: false;

		string f_niv <- bench_dossier + bench_prefixe + "_niveaux.csv";
		save ["niveau", "borne_min_um", "borne_max_um", "mineral_vise", "mineral_obtenu", "mineral_attendu",
			"mo_visee", "mo_obtenue", "mo_attendue", "pores_obtenus", "n_mineral", "n_mo", "n_pore", "n_decomposable"]
			to: f_niv format: "csv" rewrite: true header: false;
		loop r over: g.rapport {
			int i <- int(r[0]);
			save [i, r[1] / #mm * 1000, r[2] / #mm * 1000, r[3], r[5], attendus[i][0], r[4], r[6], attendus[i][1], r[7],
				int(r[8]), int(r[9]), int(r[10]), int(r[11])]
				to: f_niv format: "csv" rewrite: false header: false;
		}
	}
}

species generateur skills: [camisole] {
	list<list<float>> rapport;   // une ligne par niveau (voir l'action soil_factory_report du plugin)
	float rmse;                  // écart du tirage retenu, fraction du cube
	float duree_ms;

	init {
		location <- world.location;
		// [Initialisation] masses en g/kg (le plugin multiplie par 1000), bornes en m (le plugin convertit en cm)
		loop c over: classes_minerales {
			do granulometry(mineral_matter: c[0] #gram / #kg, min: c[1], max: c[2]);
		}
		loop c over: classes_organiques {
			do granulometry(organic_matter: c[0] #gram / #kg, min: c[1], max: c[2]);
		}

		float t0 <- gama.machine_time;
		do initialize_soil(size: soil_size, dividing_factor: division, default_species: cellule,
			number_of_try: nb_tirages, bulk_density: densite_apparente);
		duree_ms <- gama.machine_time - t0;

		rapport <- list<list<float>>(soil_factory_report());
		rmse <- float(soil_factory_rmse());
		do materialize(depth: profondeur);

		ask world { do analyser(myself); }
	}
}

species cellule skills: [apsf_particle] {
	aspect base {
		if (abs(size - soil_size / division ^ profondeur) < 1e-5 * size or my_type in ["mineral", "organic", "pore"]) {
			point p <- location + {0, 0, -size / 2};
			switch my_type {
				match "mineral" { draw cube(size) at: p color: #black border: #white; }
				match "organic" { draw cube(size) at: p color: rgb(0, 160, 60) border: #white; }
				match "pore" { draw cube(size) at: p wireframe: true color: #darkblue; }
				default { draw cube(size) at: p color: rgb(128, 128, 128, 120) border: #white; }
			}
		}
	}
}

experiment soil_factory_gui type: gui {
	parameter "Tirages" var: nb_tirages category: "SoilFactory";
	parameter "Profondeur" var: profondeur min: 1 max: 2 category: "SoilFactory";
	output {
		display vue type: opengl background: #white {
			species cellule aspect: base;
		}
	}
}

// [Bench] Batch : 5 sols pour 2 000 et 20 000 tirages. Runs séquentiels (variables static dans le plugin).
experiment soil_factory_batch type: batch repeat: 5 until: cycle >= 1 parallel: false {
	parameter "Tirages" var: nb_tirages among: [2000, 20000];
	method exhaustive;
}
