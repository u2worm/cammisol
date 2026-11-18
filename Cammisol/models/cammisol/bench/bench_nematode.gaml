/**
 * Name: cammisol_bench_nematode
 * Banc d'essai de l'effet des nématodes : même sol (grid4), plusieurs intensités de prédation.
 * Réutilise l'export de bench.gaml ; les fichiers vont dans Cammisol/bench_results/nematodes/,
 * avec le scénario dans leur nom (ex. grid4_n2_nem-B_s<graine>_global.csv).
 *
 * Les 500 agents nématodes sont gardés ; chaque agent représente plusieurs nématodes, donc seul
 * nematode_predation_rate (C prédaté par agent et par jour) change selon le scénario.
 * Prédation totale = 500 x taux ; biomasse microbienne initiale du cube ≈ 1,16 g C.
 *   actuel : 1e-7 g/j     -> 0,00005 g/j (≈ 0,004 % de la biomasse initiale par jour)
 *   A      : 1.9e-5 g/j   -> 0,0094 g/j  (≈ 0,8 %/j ; ordre de 10 nématodes par g de sol, à vérifier)
 *   B      : 2.3e-4 g/j   -> 0,116 g/j   (≈ 10 %/j)
 *   C      : 1.2e-3 g/j   -> 0,58 g/j    (≈ 50 %/j)
 */

model cammisol_bench_nematode

import "bench.gaml"

global {
	string scenario_nematodes <- "A";
	map<string, float> taux_nematodes <- [
		"actuel"::1e-7 #gram/#day,
		"A"::1.9e-5 #gram/#day,
		"B"::2.3e-4 #gram/#day,
		"C"::1.2e-3 #gram/#day
	];

	init {
		nematode_predation_rate <- taux_nematodes[scenario_nematodes];
		bench_cycles <- int(bench_jours * #day / local_step);   // recalculé au cas où la durée vient d'un paramètre
		bench_dossier <- "../../../bench_results/nematodes/";
		bench_etiquette <- "nem-" + scenario_nematodes;
		write "[Bench nématodes] Scénario " + scenario_nematodes + " : " + (nematode_predation_rate / (#gram/#day))
			+ " g C/j par agent, " + nematodes_count + " agents, soit "
			+ (nematodes_count * nematode_predation_rate / (#gram/#day)) + " g C/j au total.";
	}
}

// En interface : un scénario à la fois, avec les graphiques habituels.
experiment bench_nematode_gui parent: base_cammisol_output {
	parameter "Scénario nématodes" var: scenario_nematodes category: "Bench" among: ["actuel", "A", "B", "C"];
	parameter "Canevas" var: soil_file category: "Bench" among: ["../../../canvases/grid4.json"];
	parameter "Durée (jours)" var: bench_jours category: "Bench";
}

// Batch : 4 scénarios x 3 graines = 12 runs de 20 jours sur grid4, lancés un par un.
experiment bench_nematode type: batch repeat: 3 until: cycle > bench_cycles parallel: false {
	parameter "Scénario nématodes" var: scenario_nematodes among: ["actuel", "A", "B", "C"];
	parameter "Canevas" var: soil_file among: ["../../../canvases/grid4.json"];
	parameter "Durée (jours)" var: bench_jours among: [20.0];
	parameter "Optimisation initiale" var: do_opti among: [false];
	parameter "Pause à la fin" var: bench_pause_fin among: [false];
	method exhaustive;
}
