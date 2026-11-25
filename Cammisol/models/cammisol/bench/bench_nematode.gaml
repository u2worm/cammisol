/**
 * Name: cammisol_bench_nematode
 * Banc d'essai de l'effet des nématodes : même sol (grid4), plan croisé
 * 3 intensités de prédation x 3 quantités de DOM initial = 9 simulations (une seule graine,
 * identique pour toutes : la graine a très peu d'effet sur les résultats).
 * Réutilise l'export de bench.gaml ; fichiers dans Cammisol/bench_results/nematodes/, avec les deux
 * facteurs dans leur nom (ex. grid4_n2_pred-moyen_dom-x10_s<graine>_global.csv).
 *
 * Prédation : les 500 agents sont gardés, chaque agent représente plusieurs nématodes ; seul
 * nematode_predation_rate (C prédaté par agent et par jour) change.
 * Prédation totale = 500 x taux ; biomasse microbienne initiale du cube ≈ 1,16 g C.
 *   faible : 1.9e-5 g/j -> 0,0094 g/j (≈ 0,8 % de la biomasse initiale par jour)
 *   moyen  : 2.3e-4 g/j -> 0,116 g/j  (≈ 10 %/j)
 *   fort   : 1.2e-3 g/j -> 0,58 g/j   (≈ 50 %/j ; extinctions locales possibles)
 *
 * DOM initial : facteur appliqué ensemble au C, N et P du DOM (rapports C:N:P inchangés), par
 * rapport aux valeurs de new_grid.gaml (72 / 6 / 40 µg/g). Le DIM n'est pas modifié.
 *   x0.1 -> 7,2 µg C/g ; x1 -> 72 µg C/g ; x10 -> 720 µg C/g (≈ valeur POXC d'origine)
 */

model cammisol_bench_nematode

import "bench.gaml"

global {
	string scenario_predation <- "moyen";
	map<string, float> taux_predation <- [
		"aucune"::1e-7 #gram/#day,      // valeur d'origine du modèle (effet négligeable), pour l'interface
		"faible"::1.9e-5 #gram/#day,
		"moyen"::2.3e-4 #gram/#day,
		"fort"::1.2e-3 #gram/#day
	];
	float facteur_dom <- 10.0;

	init {
		nematode_predation_rate <- taux_predation[scenario_predation];
		bench_cycles <- int(bench_jours * #day / local_step);   // recalculé au cas où la durée vient d'un paramètre

		// Facteur sur le DOM initial. L'ordre entre cet init et celui du modèle (qui crée le sol et
		// remplit les DOM) n'est pas garanti : on modifie les concentrations (utilisées si le sol n'est
		// pas encore créé) ET les DOM déjà remplis (si le sol existe déjà). Le facteur s'applique
		// donc exactement une fois dans les deux cas.
		C_concentration_in_dom <- C_concentration_in_dom * facteur_dom;
		N_concentration_in_dom <- N_concentration_in_dom * facteur_dom;
		P_concentration_in_dom <- P_concentration_in_dom * facteur_dom;
		bool sol_deja_cree <- !empty(Dam);
		if (sol_deja_cree) {
			ask Dam { dom <- [dom[0] * facteur_dom, dom[1] * facteur_dom, dom[2] * facteur_dom]; }
		}

		bench_dossier <- "../../../bench_results/nematodes/";
		bench_etiquette <- "pred-" + scenario_predation + "_dom-x" + facteur_dom;
		write "[Bench nématodes] Prédation " + scenario_predation + " : " + (nematode_predation_rate / (#gram/#day))
			+ " g C/j par agent (" + (nematodes_count * nematode_predation_rate / (#gram/#day)) + " g C/j au total). "
			+ "DOM initial x" + facteur_dom + " : C = " + (C_concentration_in_dom / (#gram/#gram) * 1e6) + " µg/g de sol"
			+ (sol_deja_cree ? " (appliqué aux DOM déjà créés)." : " (appliqué avant la création du sol).");
	}
}

// En interface : une combinaison à la fois, avec les graphiques habituels.
experiment bench_nematode_gui parent: base_cammisol_output {
	parameter "Prédation" var: scenario_predation category: "Bench" among: ["aucune", "faible", "moyen", "fort"];
	parameter "Facteur DOM initial" var: facteur_dom category: "Bench" among: [0.1, 1.0, 10.0];
	parameter "Canevas" var: soil_file category: "Bench" among: ["../../../canvases/grid4.json"];
	parameter "Durée (jours)" var: bench_jours category: "Bench";
}

// Batch : 3 prédations x 3 DOM initiaux = 9 runs de 20 jours sur grid4, même graine, lancés un par un.
experiment bench_nematode type: batch repeat: 1 keep_seed: true until: cycle > bench_cycles parallel: false {
	parameter "Prédation" var: scenario_predation among: ["faible", "moyen", "fort"];
	parameter "Facteur DOM initial" var: facteur_dom among: [0.1, 1.0, 10.0];
	parameter "Canevas" var: soil_file among: ["../../../canvases/grid4.json"];
	parameter "Durée (jours)" var: bench_jours among: [21.0];
	parameter "Optimisation initiale" var: do_opti among: [false];
	parameter "Pause à la fin" var: bench_pause_fin among: [false];
	method exhaustive;
}
