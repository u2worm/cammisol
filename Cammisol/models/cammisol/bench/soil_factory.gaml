/**
 * Name: soil_factory_bench
 * Banc d'essai du générateur de sol SoilFactory (plugin camisole), sans biologie.
 * Génère un sol avec initialize_soil, le vérifie, le mesure niveau par niveau et l'affiche.
 *
 * Deux méthodes de génération :
 *  - "contrainte" : minéral et MO de chaque niveau égaux aux volumes visés, porosité totale = 1 - solides,
 *    pores répartis entre niveaux "au prorata des solides" ou "au hasard" ;
 *  - "sworm" : méthode d'origine du plugin (meilleur de nb_tirages tirages aléatoires).
 * Tous les tirages du plugin viennent de la graine de la simulation : même graine, même sol.
 *
 * Fichiers écrits dans Cammisol/bench_results/soil_factory/ :
 *  - generations.csv   : une ligne par run (paramètres, mesures, résultats des tests) ;
 *  - <run>_niveaux.csv : une ligne par niveau.
 * Volumes en fraction du volume du cube.
 */

model soil_factory_bench

global {
	string MINERAL <- "mineral";
	string ORGANIC <- "organic";
	string PORE <- "pore";

	// [Bench] Sol généré
	string methode <- "contrainte" among: ["contrainte", "sworm"];
	string repartition_pores <- "solides" among: ["solides", "hasard"];   // méthode "contrainte" seulement
	int nb_tirages <- 2000;            // méthode "sworm" seulement
	float soil_size <- 0.2;            // côté du cube (m)
	int division <- 10;                // cellules par côté d'un canevas
	float densite_apparente <- 1.17;   // [Initialisation] g/cm3 (ODD), passée SANS unité GAMA au plugin
	int profondeur <- 1;               // niveaux demandés à materialize
	int max_cellules <- 200000;        // au-delà, la profondeur est réduite (mémoire et affichage)
	string scenario <- "";             // batch : "contrainte_solides", "contrainte_hasard" ou "sworm"
	int copie <- 1;                    // batch de reproductibilité : deux copies d'un même run
	//float seed <- 100.0;

	// [Bench] Tolérances des tests (fractions du cube)
	float tolerance <- 1e-6;                          // calculs (volumes visés, conservation)
	float tolerance_cellule -> 0.5 / (division ^ 3);  // arrondi à une demi-cellule de canevas
	float tolerance_densite <- 0.01;                  // écart relatif de densité apparente

	// [Bench] Constantes du plugin (SoilConstant), reprises pour les calculs indépendants
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

	// [Affichage] Réglages de la vue 3D (modifiables pendant la simulation)
	int niveau_affiche <- 2;           // niveau montré : cellules de ce niveau + particules des niveaux au-dessus
	int profondeur_effective <- 1;     // niveaux réellement créés
	float coupe <- 1.0;                // fraction de la hauteur affichée (coupe horizontale)
	float opacite <- 1.0;
	bool voir_mineral <- true;
	bool voir_organic <- true;
	bool voir_agglo <- true;
	bool voir_pores <- false;
	string style_pores <- "fil de fer" among: ["fil de fer", "transparent"];
	bool voir_cadre <- true;
	float z_bas;

	// [Canevas] Contenu des canevas lu dans le plugin (canvas_cells) : canevas_cellules[n][i][j][k], niveau n
	list<string> canevas_noms <- [];
	list<list<list<list<string>>>> canevas_cellules <- [];
	list<list<int>> cellules_canevas <- [];   // pour l'affichage : [canevas, i, j, k, type]
	int canevas_affiche <- -1;                // -1 : tous les canevas côte à côte ; sinon le numéro du canevas

	// [Loupe] Fenêtre cubique du sol décrite à partir des canevas, sans agents : on descend de niveau en
	// niveau tant que le nombre de cellules à dessiner reste sous loupe_max_cellules.
	float loupe_cote <- 0.01;                 // côté de la fenêtre (m)
	float loupe_x <- 0.0;                     // position de la fenêtre dans le cube, de 0 à 1 (0 : coin du bas)
	float loupe_y <- 0.0;
	float loupe_z <- 0.0;
	int loupe_max_cellules <- 50000;
	int loupe_niveau_max <- 10;
	list<list<float>> cellules_loupe <- [];   // [x, y, z, dx, dy, dz, type], ramenées à la taille de l'affichage

	// [Bench] Sorties
	string bench_dossier <- "../../../bench_results/soil_factory/";
	string bench_etiquette <- "";      // ajoutée au nom des fichiers, vide par défaut
	string bench_graine;
	string bench_prefixe;
	string bench_cle;                  // identifie un sol : même clé et même graine -> même sol attendu

	geometry shape <- cube(soil_size);

	init {
		// [Bench] scénario du batch
		switch scenario {
			match "contrainte_solides" { methode <- "contrainte"; repartition_pores <- "solides"; }
			match "contrainte_hasard" { methode <- "contrainte"; repartition_pores <- "hasard"; }
			match "sworm" { methode <- "sworm"; }
		}
		bench_graine <- replace(string(seed), ".", "");
		string reglage <- methode = "contrainte" ? repartition_pores : "t" + nb_tirages;
		bench_cle <- methode + "_" + reglage + "_c" + (soil_size / #mm) + "mm_d" + division + "_n" + profondeur;
		bench_prefixe <- "sf_" + bench_cle + (bench_etiquette = "" ? "" : "_" + bench_etiquette) + "_s" + bench_graine;
		z_bas <- location.z - soil_size / 2;
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

	// [Bench] Plus grande profondeur <= profondeur dont le nombre de cellules reste sous max_cellules.
	// Niveau d : division^3 x (produit des cellules décomposables des canevas des niveaux 0 à d - 2).
	int profondeur_possible(list<list<float>> rapport) {
		float nb_niveau <- float(division ^ 3);
		float total <- nb_niveau;
		int d <- 1;
		loop while: d < profondeur and d < length(rapport) {
			nb_niveau <- nb_niveau * rapport[d - 1][11];
			if (total + nb_niveau > max_cellules) { break; }
			total <- total + nb_niveau;
			d <- d + 1;
		}
		return d;
	}

	// [Bench] Empreinte d'un sol : contenu des canevas + position et type de chaque cellule créée.
	// Deux runs de même clé et de même graine doivent avoir la même empreinte.
	string empreinte(list<list<float>> rapport) {
		string s <- "";
		loop r over: rapport {
			s <- s + int(r[8]) + "-" + int(r[9]) + "-" + int(r[10]) + "-" + int(r[11]) + "/";
		}
		// code du type (minéral 1, MO 2, pore 3, décomposable 4) x niveau x position dans le canevas parent
		float somme <- sum(cellule collect (
			(each.my_type = MINERAL ? 1 : (each.my_type = ORGANIC ? 2 : (each.my_type = PORE ? 3 : 4)))
			* each.niveau * (1 + each.i + division * each.j + division * division * each.k)));
		return s + string(somme);
	}

	// [Bench] Dernière empreinte enregistrée pour la même clé et la même graine ("" si aucune)
	string empreinte_precedente(string f_runs, string cle_graine) {
		string res <- "";
		if (file_exists(f_runs)) {
			list<string> lignes <- list<string>(text_file(f_runs).contents);
			if (length(lignes) > 1) {
				list<string> entete <- lignes[0] split_with ",";
				int i_cle <- entete index_of "cle_graine";
				int i_emp <- entete index_of "empreinte";
				if (i_cle < 0 or i_emp < 0) { return ""; }   // fichier d'une ancienne version
				loop n from: 1 to: length(lignes) - 1 {
					list<string> col <- lignes[n] split_with ",";
					if (length(col) > max([i_cle, i_emp]) and col[i_cle] = cle_graine) { res <- col[i_emp]; }
				}
			}
		}
		return res;
	}

	// [Affichage] Type d'une cellule : 1 minéral, 2 MO, 3 pore, 4 décomposable
	int code_type(string t) {
		return t = MINERAL ? 1 : (t = ORGANIC ? 2 : (t = PORE ? 3 : 4));
	}

	bool type_visible(int code) {
		return (code = 1 and voir_mineral) or (code = 2 and voir_organic) or (code = 3 and voir_pores) or (code = 4 and voir_agglo);
	}

	// [Affichage] Couleurs : minéral noir, MO verte, pore bleu, décomposable gris
	rgb couleur_type(int code) {
		int a <- int(opacite * 255);
		return code = 1 ? rgb(0, 0, 0, a) : (code = 2 ? rgb(0, 160, 60, a) : (code = 3 ? rgb(80, 120, 255, int(opacite * 60)) : rgb(128, 128, 128, a)));
	}

	// [Canevas] Disposition de la vue des canevas : tous sur une grille de 3 colonnes, ou un seul en grand
	float cote_canevas -> canevas_affiche < 0 ? soil_size / 3.4 : soil_size;
	point origine_canevas(int n) {
		if (canevas_affiche >= 0) { return {0.0, 0.0, z_bas}; }
		return {(n mod 3) * 1.2 * cote_canevas, (n div 3) * 1.2 * cote_canevas, z_bas};
	}

	action preparer_canevas() {
		cellules_canevas <- [];
		loop n from: 0 to: length(canevas_cellules) - 1 {
			list<list<list<string>>> cells <- canevas_cellules[n];
			loop i from: 0 to: length(cells) - 1 {
				loop j from: 0 to: length(cells[i]) - 1 {
					loop k from: 0 to: length(cells[i][j]) - 1 {
						cellules_canevas <+ [n, i, j, k, code_type(cells[i][j][k])];
					}
				}
			}
		}
	}

	// [Loupe] Indices [i0, i1, j0, j1, k0, k1] des cellules du canevas f = [canevas, x, y, z, côté] qui touchent
	// la fenêtre d'origine o et de côté w (coordonnées en m dans le cube)
	list<int> plage(list<float> f, point o, float w) {
		float c <- f[4] / division;
		return [max([0, int(floor((o.x - f[1]) / c))]), min([division - 1, int(ceil((o.x + w - f[1]) / c)) - 1]),
			max([0, int(floor((o.y - f[2]) / c))]), min([division - 1, int(ceil((o.y + w - f[2]) / c)) - 1]),
			max([0, int(floor((o.z - f[3]) / c))]), min([division - 1, int(ceil((o.z + w - f[3]) / c)) - 1])];
	}

	int nb_dans_fenetre(list<int> r) {
		return (r[1] < r[0] or r[3] < r[2] or r[5] < r[4]) ? 0 : (r[1] - r[0] + 1) * (r[3] - r[2] + 1) * (r[5] - r[4] + 1);
	}

	// [Loupe] Ajoute la partie d'une cellule (coin x, y, z, côté c) comprise dans la fenêtre, agrandie à la taille de l'affichage
	action ajouter_loupe(float x, float y, float z, float c, int code, point o, float w) {
		float x0 <- max([x, o.x]);
		float x1 <- min([x + c, o.x + w]);
		float y0 <- max([y, o.y]);
		float y1 <- min([y + c, o.y + w]);
		float z0 <- max([z, o.z]);
		float z1 <- min([z + c, o.z + w]);
		if (x1 - x0 > 1e-9 * w and y1 - y0 > 1e-9 * w and z1 - z0 > 1e-9 * w) {
			float f <- soil_size / w;
			cellules_loupe <+ [(x0 - o.x) * f, (y0 - o.y) * f, (z0 - o.z) * f, (x1 - x0) * f, (y1 - y0) * f, (z1 - z0) * f, float(code)];
		}
	}

	// [Loupe] Descente niveau par niveau dans la fenêtre. Avant chaque niveau, on compte les cellules qu'il
	// ajouterait : si le total dépasse loupe_max_cellules (ou si loupe_niveau_max est atteint), on s'arrête et
	// les cellules décomposables restantes sont dessinées en gris.
	action calculer_loupe() {
		float w <- min([loupe_cote, soil_size]);
		point o <- {loupe_x * (soil_size - w), loupe_y * (soil_size - w), loupe_z * (soil_size - w)};
		cellules_loupe <- [];
		list<list<float>> frontiere <- [[0.0, 0.0, 0.0, 0.0, soil_size]];   // cellules à développer : [canevas, x, y, z, côté]
		int niveau <- 0;
		string arret <- "plus de cellule décomposable";
		loop while: !empty(frontiere) {
			int a_ajouter <- 0;
			loop f over: frontiere { a_ajouter <- a_ajouter + nb_dans_fenetre(plage(f, o, w)); }
			if (niveau >= loupe_niveau_max or length(cellules_loupe) + a_ajouter > loupe_max_cellules) {
				arret <- niveau >= loupe_niveau_max ? "niveau max. atteint" : "budget de cellules atteint";
				loop f over: frontiere { do ajouter_loupe(x: f[1], y: f[2], z: f[3], c: f[4], code: 4, o: o, w: w); }
				break;
			}
			niveau <- niveau + 1;
			list<list<float>> suivante <- [];
			loop f over: frontiere {
				list<int> r <- plage(f, o, w);
				if (nb_dans_fenetre(r) > 0) {
					list<list<list<string>>> cells <- canevas_cellules[int(f[0])];
					float c <- f[4] / division;
					loop i from: r[0] to: r[1] {
						loop j from: r[2] to: r[3] {
							loop k from: r[4] to: r[5] {
								string t <- cells[i][j][k];
								int code <- code_type(t);
								float x <- f[1] + i * c;
								float y <- f[2] + j * c;
								float z <- f[3] + k * c;
								if (code < 4) {
									do ajouter_loupe(x: x, y: y, z: z, c: c, code: code, o: o, w: w);
								} else {
									suivante <+ [float(canevas_noms index_of t), x, y, z, c];
								}
							}
						}
					}
				}
			}
			frontiere <- suivante;
		}
		write "[Loupe] fenêtre de " + (w / #mm) + " mm au point " + {o.x / #mm, o.y / #mm, o.z / #mm} + " (mm) : "
			+ length(cellules_loupe) + " cellules, " + niveau + " niveau(x) décrit(s), plus petite cellule "
			+ (soil_size / (division ^ niveau) / #mm * 1000) + " µm (" + arret + ")";
	}

	action analyser(generateur g) {
		float cube_cm <- soil_size / #cm;
		float vol_cube <- cube_cm ^ 3;                                           // cm3
		float masse_sol <- vol_cube * densite_apparente;                         // g
		float f_min <- masse_sol / 1000 / densite_minerale / vol_cube;           // (g/kg) -> fraction du cube
		float f_mo <- masse_sol / 1000 * facteur_mo / densite_mo / vol_cube;     // (g C/kg) -> fraction du cube
		int nb <- length(g.rapport);
		bool contrainte <- methode = "contrainte";

		// [Bench] Test 1 : nombre de niveaux
		int nb_attendu <- niveaux_attendus(cube_cm);
		bool test_niveaux <- nb = nb_attendu;

		// [Bench] Test 2 : volumes visés par niveau = volumes recalculés ici
		// [Bench] Test 3 : somme des volumes visés = volume des classes comprises dans le cube
		float ecart_vise <- 0.0;
		list<list<float>> attendus <- [];
		loop r over: g.rapport {
			float lo <- r[1] / #cm;
			float hi <- min([r[2] / #cm, cube_cm]);
			float min_att <- part_classes(classes_minerales, lo, hi, cube_cm) * f_min;
			float mo_att <- part_classes(classes_organiques, lo, hi, cube_cm) * f_mo;
			attendus <+ [min_att, mo_att];
			ecart_vise <- max([ecart_vise, abs(r[3] - min_att), abs(r[4] - mo_att)]);
		}
		bool test_vise <- ecart_vise < tolerance;
		float min_vise_total <- sum(g.rapport collect each[3]);
		float mo_visee_total <- sum(g.rapport collect each[4]);
		float min_dans_cube <- part_classes(classes_minerales, 0.0, cube_cm, cube_cm) * f_min;
		float mo_dans_cube <- part_classes(classes_organiques, 0.0, cube_cm, cube_cm) * f_mo;
		float hors_cube <- sum(classes_minerales collect each[0]) * f_min - min_dans_cube
			+ sum(classes_organiques collect each[0]) * f_mo - mo_dans_cube;   // classes plus grosses que le cube
		float ecart_conservation <- max([abs(min_vise_total - min_dans_cube), abs(mo_visee_total - mo_dans_cube)]);
		bool test_conservation <- ecart_conservation < tolerance;

		// [Bench] Test 4 : cellules du niveau 1 créées par materialize = contenu du canevas racine
		list<cellule> niveau_1 <- cellule where (each.niveau = 1);
		list<int> compte_mat <- [
			length(niveau_1 where (each.my_type = MINERAL)),
			length(niveau_1 where (each.my_type = ORGANIC)),
			length(niveau_1 where (each.my_type = PORE)),
			length(niveau_1 where !(each.est_feuille))
		];
		list<int> compte_canevas <- [int(g.rapport[0][8]), int(g.rapport[0][9]), int(g.rapport[0][10]), int(g.rapport[0][11])];
		bool test_materialize <- compte_mat = compte_canevas;

		// [Bench] Mesures : porosité, volume non résolu, densité apparente du sol généré
		float porosite_visee <- 1 - min_vise_total - mo_visee_total;
		float pores_obtenus <- sum(g.rapport collect each[7]);
		float reste <- 1 - sum(g.rapport collect (each[5] + each[6] + each[7]));
		float densite_visee <- min_vise_total * densite_minerale + mo_visee_total * densite_mo;
		float densite_obtenue <- sum(g.rapport collect each[5]) * densite_minerale + sum(g.rapport collect each[6]) * densite_mo;

		// [Bench] Tests 5 à 7 (méthode "contrainte") : volumes obtenus = visés à une demi-cellule près par niveau,
		// densité apparente obtenue = visée, aucune cellule décomposable au dernier niveau
		float ecart_obtenu <- max(g.rapport collect max([abs(each[5] - each[3]), abs(each[6] - each[4])]));
		string test_volumes <- !contrainte ? "-" : (ecart_obtenu <= tolerance_cellule + tolerance ? "OK" : "ÉCHEC");
		string test_densite <- !contrainte ? "-" : (abs(densite_obtenue - densite_visee) / densite_visee < tolerance_densite ? "OK" : "ÉCHEC");
		string test_dernier <- !contrainte ? "-" : (int(last(g.rapport)[11]) = 0 ? "OK" : "ÉCHEC");

		// [Bench] Test 8 : même clé et même graine qu'un run déjà enregistré -> même empreinte
		string f_runs <- bench_dossier + "generations.csv";
		string cle_graine <- bench_cle + "_s" + bench_graine;
		string emp <- empreinte(g.rapport);
		string emp_avant <- empreinte_precedente(f_runs, cle_graine);
		string test_repro <- emp_avant = "" ? "-" : (emp_avant = emp ? "OK" : "ÉCHEC");

		write "========== SoilFactory : " + bench_prefixe + " ==========";
		write "Génération : " + (g.duree_ms / 1000) + " s ; niveaux créés : " + profondeur_effective + " (demandé " + profondeur + ")";
		write "[Test 1] niveaux : " + nb + " (attendu " + nb_attendu + ") -> " + (test_niveaux ? "OK" : "ÉCHEC");
		write "[Test 2] volumes visés = recalculés : écart max " + ecart_vise + " -> " + (test_vise ? "OK" : "ÉCHEC");
		write "[Test 3] conservation des classes dans le cube : écart " + ecart_conservation + " -> "
			+ (test_conservation ? "OK" : "ÉCHEC") + " (hors cube : " + ((100 * hors_cube) with_precision 3) + " % du cube)";
		write "[Test 4] materialize = canevas racine (minéral, MO, pore, décomposable) : " + compte_mat + " / "
			+ compte_canevas + " -> " + (test_materialize ? "OK" : "ÉCHEC");
		write "[Test 5] volumes obtenus = visés (écart max " + ecart_obtenu + ", tolérance " + tolerance_cellule + ") -> " + test_volumes;
		write "[Test 6] densité apparente obtenue " + (densite_obtenue with_precision 3) + " / visée "
			+ (densite_visee with_precision 3) + " g/cm3 -> " + test_densite;
		write "[Test 7] aucune cellule décomposable au dernier niveau -> " + test_dernier;
		write "[Test 8] reproductibilité (même clé, même graine) -> "
			+ (test_repro = "-" ? "pas de run précédent avec cette graine" : test_repro);
		// [Bench] Test 9 : contenu des canevas lu par canvas_cells = nombres de cellules du rapport du plugin
		bool test_canevas <- length(canevas_cellules) = nb;
		loop n from: 0 to: length(canevas_cellules) - 1 {
			list<int> compte <- [0, 0, 0, 0];
			loop plan over: canevas_cellules[n] {
				loop ligne over: plan {
					loop t over: ligne {
						int c <- code_type(t) - 1;
						compte[c] <- compte[c] + 1;
					}
				}
			}
			if (compte != [int(g.rapport[n][8]), int(g.rapport[n][9]), int(g.rapport[n][10]), int(g.rapport[n][11])]) {
				test_canevas <- false;
			}
		}
		write "[Test 9] canevas lus (canvas_cells) = rapport du plugin -> " + (test_canevas ? "OK" : "ÉCHEC");
		write "Niveau | cellule | minéral visé / obtenu | MO visée / obtenue | pores obtenus | canevas (min, MO, pore, déc.)";
		loop r over: g.rapport {
			int n <- int(r[0]);
			write "C" + n + " | " + (r[1] = 0.0 ? "<= " + (r[2] / #mm * 1000) : string(r[1] / #mm * 1000)) + " µm | "
				+ ((100 * r[3]) with_precision 3) + " / " + ((100 * r[5]) with_precision 3) + " % | "
				+ ((100 * r[4]) with_precision 3) + " / " + ((100 * r[6]) with_precision 3) + " % | "
				+ ((100 * r[7]) with_precision 3) + " % | "
				+ [int(r[8]), int(r[9]), int(r[10]), int(r[11])];
		}
		write "Porosité visée " + ((100 * porosite_visee) with_precision 2) + " %, pores obtenus "
			+ ((100 * pores_obtenus) with_precision 2) + " %, non résolu sous le dernier niveau "
			+ ((100 * reste) with_precision 3) + " %, RMS " + ((100 * g.rmse) with_precision 3) + " % du cube";
		write "Densité apparente de l'ODD : " + densite_apparente + " g/cm3 (la densité visée est plus forte : "
			+ "la MO s'ajoute aux 1000 g/kg de minéraux)";

		// [Bench] Export (la ligne du run est écrite après la lecture de l'empreinte précédente)
		if (!file_exists(f_runs)) {
			save ["run", "cle_graine", "graine", "methode", "repartition_pores", "nb_tirages", "cote_m", "division",
				"profondeur", "niveaux", "duree_s", "rmse", "porosite_visee", "pores_obtenus", "non_resolu",
				"densite_visee", "densite_obtenue", "empreinte", "test_niveaux", "test_vise", "test_conservation",
				"test_materialize", "test_volumes", "test_densite", "test_dernier_niveau", "test_reproductible"]
				to: f_runs format: "csv" rewrite: true header: false;
		}
		save [bench_prefixe, cle_graine, bench_graine, methode, contrainte ? repartition_pores : "-",
			contrainte ? 0 : nb_tirages, soil_size, division, profondeur_effective, nb, g.duree_ms / 1000, g.rmse,
			porosite_visee, pores_obtenus, reste, densite_visee, densite_obtenue, emp,
			test_niveaux ? "OK" : "ÉCHEC", test_vise ? "OK" : "ÉCHEC", test_conservation ? "OK" : "ÉCHEC",
			test_materialize ? "OK" : "ÉCHEC", test_volumes, test_densite, test_dernier, test_repro]
			to: f_runs format: "csv" rewrite: false header: false;

		string f_niv <- bench_dossier + bench_prefixe + "_niveaux.csv";
		save ["niveau", "borne_min_um", "borne_max_um", "mineral_vise", "mineral_obtenu", "mineral_attendu",
			"mo_visee", "mo_obtenue", "mo_attendue", "pores_obtenus", "n_mineral", "n_mo", "n_pore", "n_decomposable"]
			to: f_niv format: "csv" rewrite: true header: false;
		loop r over: g.rapport {
			int n <- int(r[0]);
			save [n, r[1] / #mm * 1000, r[2] / #mm * 1000, r[3], r[5], attendus[n][0], r[4], r[6], attendus[n][1], r[7],
				int(r[8]), int(r[9]), int(r[10]), int(r[11])]
				to: f_niv format: "csv" rewrite: false header: false;
		}
	}
}

species generateur skills: [camisole] {
	list<list<float>> rapport;   // une ligne par niveau (voir l'action soil_factory_report du plugin)
	float rmse;                  // écart entre volumes visés et obtenus, fraction du cube
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

		// [Initialisation] noms des options du plugin
		string generation <- methode = "contrainte" ? "constrained" : "random_search";
		string pores <- repartition_pores = "solides" ? "solid" : "random";

		float t0 <- gama.machine_time;
		do initialize_soil(size: soil_size, dividing_factor: division, default_species: cellule,
			number_of_try: nb_tirages, bulk_density: densite_apparente, generation: generation, pore_distribution: pores);
		duree_ms <- gama.machine_time - t0;

		rapport <- list<list<float>>(soil_factory_report());
		rmse <- float(soil_factory_rmse());

		// [Canevas] lecture du contenu de chaque canevas (C0 = racine)
		loop n from: 0 to: length(rapport) - 1 {
			canevas_noms <+ "C" + n;
			canevas_cellules <+ list<list<list<string>>>(canvas_cells(name: "C" + n));
		}

		ask world { profondeur_effective <- profondeur_possible(myself.rapport); }
		if (profondeur_effective < profondeur) {
			warn "[Bench] profondeur réduite à " + profondeur_effective + " : plus de " + max_cellules + " cellules sinon";
		}
		do materialize(depth: profondeur_effective);

		ask world {
			do analyser(myself);
			do preparer_canevas();
			do calculer_loupe();
		}
	}

	aspect cadre {
		if (voir_cadre) { draw cube(soil_size) at: location + {0, 0, -soil_size / 2} wireframe: true color: #red; }
	}
}

species cellule skills: [apsf_particle] {
	// [Affichage] niveau de la cellule (1 = cellules du canevas racine), déduit de sa taille
	int niveau -> int(round(ln(soil_size / size) / ln(division)));
	bool est_feuille -> my_type in [MINERAL, ORGANIC, PORE];

	// [Affichage] cellules du niveau montré, plus les particules (feuilles) des niveaux plus grossiers
	aspect base {
		int vu <- min([niveau_affiche, profondeur_effective]);
		if ((niveau = vu or (niveau < vu and est_feuille)) and location.z <= z_bas + coupe * soil_size) {
			point p <- location + {0, 0, -size / 2};
			int alpha <- int(opacite * 255);
			switch my_type {
				match MINERAL {
					if (voir_mineral) { draw cube(size) at: p color: rgb(0, 0, 0, alpha) border: #white; }
				}
				match ORGANIC {
					if (voir_organic) { draw cube(size) at: p color: rgb(0, 160, 60, alpha) border: #white; }
				}
				match PORE {
					if (voir_pores) {
						if (style_pores = "fil de fer") {
							draw cube(size) at: p wireframe: true color: #darkblue;
						} else {
							draw cube(size) at: p color: rgb(80, 120, 255, int(opacite * 60)) border: #darkblue;
						}
					}
				}
				default {
					if (voir_agglo) { draw cube(size) at: p color: rgb(128, 128, 128, alpha) border: #white; }
				}
			}
		}
	}
}

experiment soil_factory_gui type: gui {
	parameter "Graine" var: seed init: 100.0 category: "Sol";
	parameter "Méthode" var: methode category: "Sol";
	parameter "Répartition des pores (contrainte)" var: repartition_pores category: "Sol";
	parameter "Tirages (sworm)" var: nb_tirages category: "Sol";
	parameter "Côté du cube (m)" var: soil_size among: [0.2, 0.02, 0.002] category: "Sol";
	parameter "Cellules par côté" var: division among: [10, 5] category: "Sol";
	parameter "Profondeur" var: profondeur min: 1 max: 3 category: "Sol";
	parameter "Cellules max." var: max_cellules category: "Sol";

	parameter "Niveau affiché" var: niveau_affiche min: 1 max: 3 category: "Affichage";
	parameter "Coupe (fraction de hauteur)" var: coupe min: 0.05 max: 1.0 step: 0.05 category: "Affichage";
	parameter "Opacité" var: opacite min: 0.05 max: 1.0 step: 0.05 category: "Affichage";
	parameter "Minéral (noir)" var: voir_mineral category: "Affichage";
	parameter "Matière organique (vert)" var: voir_organic category: "Affichage";
	parameter "Décomposable (gris)" var: voir_agglo category: "Affichage";
	parameter "Pores (bleu)" var: voir_pores category: "Affichage";
	parameter "Style des pores" var: style_pores category: "Affichage";
	parameter "Cadre (rouge)" var: voir_cadre category: "Affichage";

	parameter "Canevas affiché (-1 : tous)" var: canevas_affiche min: -1 max: 9 category: "Canevas";

	parameter "Côté de la fenêtre (m)" var: loupe_cote among: [0.2, 0.02, 0.01, 0.002, 0.001, 0.0002, 0.0001] category: "Loupe";
	parameter "Position x (0 à 1)" var: loupe_x min: 0.0 max: 1.0 step: 0.05 category: "Loupe";
	parameter "Position y (0 à 1)" var: loupe_y min: 0.0 max: 1.0 step: 0.05 category: "Loupe";
	parameter "Position z (0 : bas)" var: loupe_z min: 0.0 max: 1.0 step: 0.05 category: "Loupe";
	parameter "Cellules max. (loupe)" var: loupe_max_cellules category: "Loupe";
	parameter "Niveau max. (loupe)" var: loupe_niveau_max min: 1 max: 12 category: "Loupe";
	// [Loupe] la fenêtre est recalculée à la demande, après un changement des paramètres de la loupe
	user_command "Recalculer la loupe" { ask simulation { do calculer_loupe(); } }

	output {
		// [Affichage] Sol créé par materialize (agents)
		//display vue type: opengl background: #white {
		//	species generateur aspect: cadre;
		//	species cellule aspect: base;
		//}

		// [Canevas] Chaque canevas en 3D : une cellule décomposable d'un niveau est une copie du canevas suivant
		display canevas type: opengl background: #white {
			graphics "canevas" {
				float s <- cote_canevas;
				float c <- s / division;
				loop e over: cellules_canevas {
					if ((canevas_affiche < 0 or e[0] = canevas_affiche) and type_visible(e[4]) and e[3] * c <= coupe * s) {
						point o <- origine_canevas(e[0]);
						bool fil <- e[4] = 3 and style_pores = "fil de fer";
						draw cube(c) at: o + {(e[1] + 0.5) * c, (e[2] + 0.5) * c, e[3] * c}
							wireframe: fil color: fil ? #darkblue : couleur_type(e[4]) border: e[4] = 3 ? #darkblue : #white;
					}
				}
				loop n from: 0 to: length(canevas_noms) - 1 {
					if (canevas_affiche < 0 or n = canevas_affiche) {
						point o <- origine_canevas(n);
						if (voir_cadre) { draw cube(s) at: o + {s / 2, s / 2, 0} wireframe: true color: #red; }
						draw canevas_noms[n] + " : cellules de " + (soil_size / (division ^ (n + 1)) / #mm * 1000) + " µm"
							at: o + {0, -0.08 * s, s} color: #black;
					}
				}
			}
		}

		// [Loupe] Fenêtre du sol agrandie à la taille de l'affichage (cadre rouge = bords de la fenêtre)
		display loupe type: opengl background: #white {
			graphics "loupe" {
				loop e over: cellules_loupe {
					int code <- int(e[6]);
					if (type_visible(code) and e[2] <= coupe * soil_size) {
						bool fil <- code = 3 and style_pores = "fil de fer";
						draw box(e[3], e[4], e[5]) at: {e[0] + e[3] / 2, e[1] + e[4] / 2, z_bas + e[2]}
							wireframe: fil color: fil ? #darkblue : couleur_type(code) border: code = 3 ? #darkblue : #white;
					}
				}
				if (voir_cadre) { draw cube(soil_size) at: {soil_size / 2, soil_size / 2, z_bas} wireframe: true color: #red; }
			}
		}
	}
}

// [Bench] Comparaison des méthodes : 5 sols par scénario (graines différentes). Runs séquentiels :
// le plugin partage des variables static (flux aléatoire, taille des canevas).
experiment soil_factory_batch type: batch repeat: 5 until: cycle >= 1 parallel: false {
	parameter "Scénario" var: scenario among: ["contrainte_solides", "contrainte_hasard", "sworm"];
	method exhaustive;
}

// [Bench] Reproductibilité : pour chaque scénario, deux copies avec la même graine (keep_seed).
// La seconde doit afficher [Test 8] OK.
experiment soil_factory_repro type: batch repeat: 1 keep_seed: true until: cycle >= 1 parallel: false {
	parameter "Scénario" var: scenario among: ["contrainte_solides", "contrainte_hasard", "sworm"];
	parameter "Copie" var: copie among: [1, 2];
	method exhaustive;
}
