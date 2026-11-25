/**
 * Name: cammisol_bench
 * Banc d'essai : lance CAMMISOL sur des canevas fixes et exporte les sorties en CSV.
 * Ce fichier ne modifie pas le modèle : il l'importe et ajoute seulement l'export et les expériences.
 *
 * Fichiers écrits dans Cammisol/bench_results/ (un jeu par run, nommé <canevas>_n<niveaux>_s<graine>) :
 *  - bench_runs.csv          : une ligne par run (composition du sol, paramètres) ;
 *  - <run>_global.csv        : totaux du domaine, une ligne par jour ;
 *  - <run>_pores.csv         : une ligne par jour et par taille de pore ;
 *  - <run>_organiques.csv    : une ligne par jour et par taille d'organique.
 * Masses en grammes.
 */

model cammisol_bench

import "../cammisol.gaml"

global {
	// [Bench] Durée, fréquence d'export, dossier de sortie (même profondeur que new_grid.gaml)
	float bench_jours <- 10.0;
	int bench_cycles <- int(bench_jours * #day / local_step);
	int bench_periode <- 12;                       // pas entre deux lignes sur les csv
	string bench_dossier <- "../../../bench_results/";
	bool bench_pause_fin <- true;                  // pause à la fin (interface) ; false en batch
	string bench_etiquette <- "";                  // ajoutée au nom des fichiers (ex. scénario), vide par défaut

	// [Bench] État interne
	bool bench_pret <- false;
	bool bench_apres_pas <- false;                 // l'export s'exécute-t-il après le reflex principal ?
	string bench_canevas;
	string bench_graine;
	string bench_prefixe;
	float bench_dernier_temps;
	int bench_dernier_cycle;
	list<float> bench_tailles_pores;
	list<list<PoreParticle>> bench_groupes_pores;
	list<float> bench_tailles_org;
	list<list<OrganicParticle>> bench_groupes_org;
	list<float> bench_C_init_org;                  // C POM initial par taille (répartition de init_POM)

	reflex bench_export when: cycle mod bench_periode = 0 {
		if (!bench_pret) { do bench_preparer(); }
		do bench_ecrire();
	}

	reflex bench_fin when: bench_pause_fin and cycle >= bench_cycles {
		write "[Bench] Fin : " + bench_jours + " jours simulés. Fichiers : " + bench_dossier + bench_prefixe + "_*.csv";
		do pause;
	}

	action bench_preparer() {
		// mb_ensure_reference est appelée au début du reflex principal : si la référence est déjà
		// capturée au cycle 0, l'export s'exécute après le pas (les valeurs sont celles de fin de pas).
		bench_apres_pas <- mb_reference_captured;
		bench_canevas <- replace(last(soil_file split_with "/"), ".json", "");
		bench_graine <- replace(string(seed), ".", "");
		bench_prefixe <- bench_canevas + "_n" + niveaux + (bench_etiquette = "" ? "" : "_" + bench_etiquette) + "_s" + bench_graine;

		// Classes de taille (la structure ne change pas pendant la simulation)
		list<PoreParticle> pores <- list(PoreParticle);
		bench_tailles_pores <- remove_duplicates(pores collect each.size) sort_by (- each);
		bench_groupes_pores <- [];
		loop s over: bench_tailles_pores { bench_groupes_pores <+ pores where (abs(each.size - s) <= 1e-9 * s); }
		list<OrganicParticle> orgs <- OrganicParticle where (!each.in_pore);
		bench_tailles_org <- remove_duplicates(orgs collect each.size) sort_by (- each);
		bench_groupes_org <- [];
		loop s over: bench_tailles_org { bench_groupes_org <+ orgs where (abs(each.size - s) <= 1e-9 * s); }
		float v_org <- sum(orgs collect (each.size ^ 3));
		bench_C_init_org <- [];
		loop g over: bench_groupes_org {
			bench_C_init_org <+ C_concentration_in_pom * soil_weight * sum(g collect (each.size ^ 3)) / v_org;
		}

		// Description du run
		float v_dom <- soil_size ^ 3;
		string f_runs <- bench_dossier + "bench_runs.csv";
		if (!file_exists(f_runs)) {
			save ["run", "canevas", "niveaux", "graine", "nb_pores", "nb_organiques", "nb_mineraux",
				"vol_pores", "vol_organiques", "vol_mineraux", "nb_tailles_pores", "nb_tailles_org",
				"nematodes", "do_opti", "voisinage_faces", "export_apres_pas"]
				to: f_runs format: "csv" rewrite: true header: false;
		}
		save [bench_prefixe, bench_canevas, niveaux, bench_graine,
			length(pores), length(orgs), length(MineralParticle),
			sum(pores collect (each.size ^ 3)) / v_dom, v_org / v_dom,
			sum(MineralParticle collect (each.size ^ 3)) / v_dom,
			length(bench_tailles_pores), length(bench_tailles_org),
			length(Nematode), do_opti, voisinage_faces, bench_apres_pas]
			to: f_runs format: "csv" rewrite: false header: false;

		// En-têtes des fichiers du run
		save ["jour", "pas", "ms_par_pas",
			"C_pom_labile", "C_pom_recal", "N_pom", "P_pom",
			"C_dom", "N_dom", "P_dom", "N_dim", "P_dim",
			"C_O", "C_F", "C_M", "C_cytosol", "actif_O", "actif_F", "actif_M",
			"CO2_microbes", "CO2_nematodes", "nematodes_eveilles",
			"mb_C_err_max", "mb_N_err_max", "mb_P_err_max", "nb_negatifs"]
			to: bench_dossier + bench_prefixe + "_global.csv" format: "csv" rewrite: true header: false;
		save ["jour", "taille_mm", "nb", "part_volume_pores", "C_microbes", "C_dom", "N_dim", "P_dim",
			"actif_moyen", "nematodes", "nematodes_eveilles"]
			to: bench_dossier + bench_prefixe + "_pores.csv" format: "csv" rewrite: true header: false;
		save ["jour", "taille_mm", "nb", "C_pom_initial", "C_pom", "C_pom_labile", "N_pom", "P_pom"]
			to: bench_dossier + bench_prefixe + "_organiques.csv" format: "csv" rewrite: true header: false;

		write "[Bench] Run " + bench_prefixe + " : " + length(pores) + " pores (" + length(bench_tailles_pores)
			+ " tailles), " + length(orgs) + " organiques (" + length(bench_tailles_org) + " tailles), export "
			+ (bench_apres_pas ? "après" : "avant") + " le pas principal, durée " + bench_jours + " jours.";
		bench_dernier_temps <- gama.machine_time;
		bench_dernier_cycle <- cycle;
		bench_pret <- true;
	}

	action bench_ecrire() {
		int pas <- cycle + (bench_apres_pas ? 1 : 0);
		float jour <- pas * local_step / #day;
		float ms <- cycle > (bench_dernier_cycle) ? (gama.machine_time - bench_dernier_temps) / (cycle - bench_dernier_cycle) : 0.0;
		bench_dernier_temps <- gama.machine_time;
		bench_dernier_cycle <- cycle;

		// Totaux du domaine
		float C_O <- sum(O_Strategist collect each.C);
		float C_F <- sum(F_Strategist collect each.C);
		float C_M <- sum(M_Strategist collect each.C);
		save [jour, pas, ms,
			sum(OrganicParticle collect each.C_labile) / #gram,
			sum(OrganicParticle collect each.C_recalcitrant) / #gram,
			sum(OrganicParticle collect (each.N_labile + each.N_recalcitrant)) / #gram,
			sum(OrganicParticle collect (each.P_labile + each.P_recalcitrant)) / #gram,
			sum(Dam collect each.dom[2]) / #gram, sum(Dam collect each.dom[0]) / #gram,
			sum(Dam collect each.dom[1]) / #gram, sum(Dam collect each.dim[0]) / #gram,
			sum(Dam collect each.dim[1]) / #gram,
			C_O / #gram, C_F / #gram, C_M / #gram,
			sum(mb_microbes() collect each.cytosol_C) / #gram,
			C_O > 0 ? sum(O_Strategist collect (each.C * each.active_rate)) / C_O : 0.0,
			C_F > 0 ? sum(F_Strategist collect (each.C * each.active_rate)) / C_F : 0.0,
			C_M > 0 ? sum(M_Strategist collect (each.C * each.active_rate)) / C_M : 0.0,
			microbe_CO2_emissions / #gram, nematode_CO2_emissions / #gram,
			length(Nematode where each.awake),
			mb_C_max_rel_error, mb_N_max_rel_error, mb_P_max_rel_error, mb_negative_count]
			to: bench_dossier + bench_prefixe + "_global.csv" format: "csv" rewrite: false header: false;

		// Par taille de pore
		loop i from: 0 to: length(bench_tailles_pores) - 1 {
			list<PoreParticle> g <- bench_groupes_pores[i];
			list<MicrobePopulation> pops <- g accumulate each.populations;
			float C_mic <- sum(pops collect each.C);
			list<Nematode> nems <- Nematode where (g contains each.current_pore);
			save [jour, bench_tailles_pores[i] / #mm, length(g), sum(g collect each.part_volume),
				C_mic / #gram, sum(g collect each.dam.dom[2]) / #gram,
				sum(g collect each.dam.dim[0]) / #gram, sum(g collect each.dam.dim[1]) / #gram,
				C_mic > 0 ? sum(pops collect (each.C * each.active_rate)) / C_mic : 0.0,
				length(nems), length(nems where each.awake)]
				to: bench_dossier + bench_prefixe + "_pores.csv" format: "csv" rewrite: false header: false;
		}

		// Par taille d'organique
		loop i from: 0 to: length(bench_tailles_org) - 1 {
			list<OrganicParticle> g <- bench_groupes_org[i];
			save [jour, bench_tailles_org[i] / #mm, length(g), bench_C_init_org[i] / #gram,
				sum(g collect (each.C_labile + each.C_recalcitrant)) / #gram,
				sum(g collect each.C_labile) / #gram,
				sum(g collect (each.N_labile + each.N_recalcitrant)) / #gram,
				sum(g collect (each.P_labile + each.P_recalcitrant)) / #gram]
				to: bench_dossier + bench_prefixe + "_organiques.csv" format: "csv" rewrite: false header: false;
		}
	}
}

// [Bench] En interface : un run, avec les graphiques habituels. Choisir le canevas et le nombre de niveaux.
experiment bench_gui parent: base_cammisol_output {
	parameter "Canevas" var: soil_file category: "Bench" among: [
		"../../../canvases/grid4.json", "../../../canvases/grid2.json", "../../../canvases/grid5_fractal_test.json"];
	parameter "Niveaux" var: niveaux category: "Bench" min: 1 max: 5;
	parameter "Durée (jours)" var: bench_jours category: "Bench";
}

// [Bench] Batch : grid5_fractal_test à 1, 2 et 3 niveaux, 3 graines chacun (9 runs).
// Un seul canevas par batch et runs séquentiels : Apsf.DIVISION_PER_LEVEL est static dans le plugin.
experiment bench_grid5 type: batch repeat: 3 until: cycle > bench_cycles parallel: false {
	parameter "Canevas" var: soil_file among: ["../../../canvases/grid5_fractal_test.json"];
	parameter "Niveaux" var: niveaux among: [1, 2, 3];
	parameter "Optimisation initiale" var: do_opti among: [false];
	parameter "Pause à la fin" var: bench_pause_fin among: [false];
	method exhaustive;
}
