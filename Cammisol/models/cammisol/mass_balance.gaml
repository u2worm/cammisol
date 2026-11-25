/**
 * Name: mass_balance
 * Auteur : T. Bonetti
 *
 * Contrôle du bilan de masse C / N / P à chaque pas de temps, et test de positivité
 * des stocks.
 *
 * PRINCIPE
 * --------
 * Sauf entrées ou sorties explicites, le modèle ne doit ni créer ni détruire de matière.
 * Pour chaque élément X dans {C, N, P}, on vérifie à chaque pas :
 *
 *     Erreur_X(t) = [Stocks_X(t) + Sorties_X cumulées(t)] - [Stocks_X(0) + Entrées_X cumulées(t)]
 *     Erreur relative_X(t) = Erreur_X(t) / (Stocks_X(0) + Entrées_X cumulées(t))
 *
 *
 * COMPARTIMENTS COMPTÉS (stocks)
 * ------------------------------
 *   - POM    : OrganicParticle (labile + récalcitrant)                 -> C, N, P
 *   - DOM    : Dam.dom  (indices : 0 = N, 1 = P, 2 = C)                -> C, N, P
 *   - DIM    : Dam.dim  (indices : 0 = N, 1 = P)                       -> N, P
 *   - Microbes : toutes les populations (O, F, M) : structural + cytosol -> C, N, P
 *   - Nématodes : contenu de l'estomac (stomack_C/N/P)                 -> C, N, P
 *     (la biomasse propre des nématodes est constante dans le modèle, elle n'est pas comptée)
 *   - MineralParticle : N, P (nuls actuellement, comptés par précaution pour la future MAM)
 *
 
 *  l'action mb_register_input() doit etre appelée AVANT d'ajouter la matière dans les compartiments.
 * UTILISATION
 * -----------
 *   - Appeler do mb_ensure_reference() au début du pas, et do mb_check() à la fin du pas
 */

model mass_balance

import "environment/new_grid.gaml"
import "nematode/nematode.gaml"

global {
	// =========================================================================================
	// Paramètres (modifiables dans l'interface, catégorie "Bilan de masse")
	// =========================================================================================

	/** Active / désactive le contrôle du bilan de masse. */
	bool mb_enabled <- true;

	/** Seuil d'alerte sur la valeur absolue de l'erreur relative (C, N et P). */
	float mb_tolerance <- 1e-9;

	/** Met la simulation en pause lorsqu'une alerte (bilan ou positivité) est levée. */
	bool mb_pause_on_error <- false;

	/** Active le test de positivité des stocks. */
	bool mb_check_positivity <- true;

	/**
	 * Tolérance du test de positivité, relative au stock initial de l'élément.
	 */
	float mb_positivity_tolerance <- 1e-12;

	/** Nombre maximal de messages d'avertissement écrits dans la console (évite de la saturer). */
	int mb_max_warnings <- 30;

	// =========================================================================================
	// État interne du contrôle (lu par les graphiques)
	// =========================================================================================

	/** Vrai une fois l'état de référence capturé (au premier pas de simulation). */
	bool mb_reference_captured <- false;

	/** Stocks de référence (+ sorties déjà cumulées au moment de la capture, normalement 0). */
	float mb_C_ref <- 0.0;
	float mb_N_ref <- 0.0;
	float mb_P_ref <- 0.0;

	/** Entrées extérieures cumulées depuis la capture de référence. */
	float mb_C_input <- 0.0;
	float mb_N_input <- 0.0;
	float mb_P_input <- 0.0;

	/** Erreurs absolues (kg) et relatives (sans unité) du dernier contrôle. Signe + = création. */
	float mb_C_error <- 0.0;
	float mb_N_error <- 0.0;
	float mb_P_error <- 0.0;
	float mb_C_rel_error <- 0.0;
	float mb_N_rel_error <- 0.0;
	float mb_P_rel_error <- 0.0;

	/** Plus grande erreur relative (en valeur absolue) observée depuis le début. */
	float mb_C_max_rel_error <- 0.0;
	float mb_N_max_rel_error <- 0.0;
	float mb_P_max_rel_error <- 0.0;

	/** Détail des stocks de C au dernier contrôle (kg), pour le graphique de répartition. */
	float mb_C_pom <- 0.0;
	float mb_C_dom <- 0.0;
	float mb_C_microbes <- 0.0;
	float mb_C_nematodes <- 0.0;
	float mb_C_co2 <- 0.0;


	/** meme chose pour le P */
	float mb_P_pom <- 0.0;
	float mb_P_dam <- 0.0;
	float mb_P_microbes <- 0.0;
	float mb_P_nematodes <- 0.0;
	float mb_P_mam <- 0.0;
	
	/** meme chose pour le N */
	float mb_N_pom <- 0.0;
	float mb_N_dam <- 0.0;
	float mb_N_microbes <- 0.0;
	float mb_N_nematodes <- 0.0;
	float mb_N_mam <- 0.0;
	
	
	/** Nombre de valeurs négatives détectées au dernier contrôle, et depuis le début. */
	int mb_negative_count <- 0;
	int mb_negative_total <- 0;

	/** Nombre de messages déjà écrits (plafonné par mb_max_warnings). */
	int mb_warning_count <- 0;

	// =========================================================================================
	// Calcul des stocks
	// =========================================================================================

	/**
	 * Toutes les populations microbiennes, quelle que soit leur strategie.
	 * On utilise of_generic_species car le nom de l'espèce mère seul (MicrobePopulation)
	 * ne renvoie pas les agents des sous-espèces.
	 */
	list<MicrobePopulation> mb_microbes() {
		return list<MicrobePopulation>(agents of_generic_species MicrobePopulation);
	}

	// --- Carbone ---
	float mb_stock_C_pom() { return sum(OrganicParticle collect (each.C_labile + each.C_recalcitrant)); }
	float mb_stock_C_dom() { return sum(Dam collect each.dom[2]); }
	float mb_stock_C_microbes() { return sum(mb_microbes() collect (each.C + each.cytosol_C)); }
	float mb_stock_C_nematodes() { return sum(Nematode collect each.stomack_C); }

	float mb_stock_C() {
		return mb_stock_C_pom() + mb_stock_C_dom() + mb_stock_C_microbes() + mb_stock_C_nematodes();
	}

	// --- Azote ---
	float mb_stock_N_pom() { return sum(OrganicParticle collect (each.N_labile + each.N_recalcitrant)) ;}
	float mb_stock_N_dam() { return sum(Dam collect (each.dom[0] + each.dim[0])) ;}
	float mb_stock_N_microbes() { return sum(mb_microbes() collect (each.N + each.cytosol_N)) ;}
	float mb_stock_N_nematodes() { return sum(Nematode collect each.stomack_N) ;}
	float mb_stock_N_mam() { return sum(MineralParticle collect each.N) ;}
	
	
	float mb_stock_N() {
		return mb_stock_N_pom() + mb_stock_N_dam() + mb_stock_N_microbes() +  mb_stock_N_nematodes() +  mb_stock_N_mam();
	}

	// --- Phosphore ---
	float mb_stock_P_pom() { return sum(OrganicParticle collect (each.P_labile + each.P_recalcitrant)) ;}
	float mb_stock_P_dam() { return sum(Dam collect (each.dom[1] + each.dim[1])) ;}
	float mb_stock_P_microbes() { return sum(mb_microbes() collect (each.P + each.cytosol_P)) ;}
	float mb_stock_P_nematodes() { return sum(Nematode collect each.stomack_P) ;}
	float mb_stock_P_mam() { return sum(MineralParticle collect each.P) ;}
	
	
	float mb_stock_P() {
		return mb_stock_P_pom() + mb_stock_P_dam() + mb_stock_P_microbes() + mb_stock_P_nematodes() + mb_stock_P_mam();
	}

	// --- Sorties cumulées ---
	/** C sorti du système : CO2 respiré par les microbes et les nématodes. */
	float mb_output_C() { return microbe_CO2_emissions + nematode_CO2_emissions; }
	/** Aucune sortie de N pour l'instant (à compléter : prélèvement des plantes, lessivage...). */
	float mb_output_N() { return 0.0; }
	/** Aucune sortie de P pour l'instant (à compléter : prélèvement des plantes, lessivage...). */
	float mb_output_P() { return 0.0; }

	// =========================================================================================
	// Actions publiques
	// =========================================================================================

	/**
	 * Capture l'état de référence si ce n'est pas déjà fait.
	 * Appelée au début du premier pas (et non dans init) car certaines expériences
	 * (ex. test_microbes) modifient le POM dans leur propre init, exécuté APRÈS celui du modèle.
	 */
	action mb_ensure_reference() {
		if (!mb_enabled or mb_reference_captured) { return; }
		mb_C_ref <- mb_stock_C() + mb_output_C();
		mb_N_ref <- mb_stock_N() + mb_output_N();
		mb_P_ref <- mb_stock_P() + mb_output_P();
		mb_reference_captured <- true;
		write "[Bilan de masse] Référence capturée (cycle " + cycle + ") : C = " + (mb_C_ref / #gram)
			+ " g, N = " + (mb_N_ref / #gram) + " g, P = " + (mb_P_ref / #gram) + " g";
	}

	/**
	 * Déclare un apport extérieur de matière (engrais, apport test...).
	 * DOIT être appelée AVANT d'ajouter la matière dans les compartiments : la référence est
	 * capturée ici si nécessaire, pour que l'apport ne soit pas compté deux fois.
	 */
	action mb_register_input (float C, float N, float P) {
		if (!mb_enabled) { return; }
		do mb_ensure_reference();
		mb_C_input <- mb_C_input + C;
		mb_N_input <- mb_N_input + N;
		mb_P_input <- mb_P_input + P;
	}

	/**
	 * Contrôle complet : bilan C/N/P puis (optionnellement) positivité.
	 * À appeler à la fin de chaque pas, une fois tous les processus exécutés.
	 */
	action mb_check() {
		if (!mb_enabled) { return; }
		do mb_ensure_reference();

		// Détail du C (stocké pour le graphique de répartition)
		mb_C_pom <- mb_stock_C_pom();
		mb_C_dom <- mb_stock_C_dom();
		mb_C_microbes <- mb_stock_C_microbes();
		mb_C_nematodes <- mb_stock_C_nematodes();
		mb_C_co2 <- mb_output_C();
		
		// meme chose pour le N
		mb_N_pom <- mb_stock_N_pom();
		mb_N_dam <- mb_stock_N_dam();
		mb_N_microbes <- mb_stock_N_microbes();
		mb_N_nematodes <- mb_stock_N_nematodes();
		mb_N_mam <- mb_stock_N_mam();
		
		// meme chose pour le P
		mb_P_pom <- mb_stock_P_pom();
		mb_P_dam <- mb_stock_P_dam();
		mb_P_microbes <- mb_stock_P_microbes();
		mb_P_nematodes <- mb_stock_P_nematodes();
		mb_P_mam <- mb_stock_P_mam();

		// Bilan : (stocks + sorties) - (référence + entrées)
		float C_expected <- mb_C_ref + mb_C_input;
		float N_expected <- mb_N_ref + mb_N_input;
		float P_expected <- mb_P_ref + mb_P_input;

		mb_C_error <- (mb_C_pom + mb_C_dom + mb_C_microbes + mb_C_nematodes + mb_C_co2) - C_expected;
		mb_N_error <- (mb_stock_N() + mb_output_N()) - N_expected;
		mb_P_error <- (mb_stock_P() + mb_output_P()) - P_expected;

		// Erreur relative (garde-fou si la quantité attendue est nulle)
		mb_C_rel_error <- C_expected != 0.0 ? mb_C_error / C_expected : 0.0;
		mb_N_rel_error <- N_expected != 0.0 ? mb_N_error / N_expected : 0.0;
		mb_P_rel_error <- P_expected != 0.0 ? mb_P_error / P_expected : 0.0;

		mb_C_max_rel_error <- max(mb_C_max_rel_error, abs(mb_C_rel_error));
		mb_N_max_rel_error <- max(mb_N_max_rel_error, abs(mb_N_rel_error));
		mb_P_max_rel_error <- max(mb_P_max_rel_error, abs(mb_P_rel_error));

		bool alert <- false;
		if (abs(mb_C_rel_error) > mb_tolerance) { alert <- true; do mb_warn("C", mb_C_rel_error, mb_C_error); }
		if (abs(mb_N_rel_error) > mb_tolerance) { alert <- true; do mb_warn("N", mb_N_rel_error, mb_N_error); }
		if (abs(mb_P_rel_error) > mb_tolerance) { alert <- true; do mb_warn("P", mb_P_rel_error, mb_P_error); }

		if (mb_check_positivity) {
			do mb_positivity();
			if (mb_negative_count > 0) { alert <- true; }
		}

		if (alert and mb_pause_on_error) {
			write "[Bilan de masse] Simulation mise en pause au cycle " + cycle + " (mb_pause_on_error = true)";
			do pause;
		}
	}

	// =========================================================================================
	// Actions internes
	// =========================================================================================

	/** Écrit un avertissement de bilan, dans la limite de mb_max_warnings messages. */
	action mb_warn (string element, float rel_error, float abs_error) {
		if (mb_warning_count < mb_max_warnings) {
			write "[Bilan de masse] ALERTE " + element + " au cycle " + cycle
				+ " : erreur relative = " + rel_error
				+ " (absolue = " + (abs_error / #gram) + " g ; "
				+ (abs_error > 0 ? "création" : "destruction") + " de matière)";
			mb_warning_count <- mb_warning_count + 1;
			if (mb_warning_count = mb_max_warnings) {
				write "[Bilan de masse] Nombre maximal de messages atteint, les suivants ne sont plus affichés.";
			}
		}
	}

	/**
	 * Test de positivité : parcourt tous les compartiments et compte les valeurs négatives
	 * (au-delà de la tolérance). Les premières sont détaillées dans la console avec le nom
	 * de l'agent et du compartiment concerné.
	 * Convention : chaque clé commence par l'élément concerné (C, N ou P), ce qui permet
	 * d'appliquer la tolérance propre à cet élément.
	 */
	action mb_positivity() {
		float tol_C <- mb_positivity_tolerance * abs(mb_C_ref);
		float tol_N <- mb_positivity_tolerance * abs(mb_N_ref);
		float tol_P <- mb_positivity_tolerance * abs(mb_P_ref);

		// Liste des (agent, compartiment, valeur) à contrôler
		list<agent> owners <- [];
		list<string> mb_keys <- [];
		list<float> mb_values <- [];

		ask OrganicParticle {
			owners <+ self; mb_keys <+ "C_labile"; mb_values <+ C_labile;
			owners <+ self; mb_keys <+ "N_labile"; mb_values <+ N_labile;
			owners <+ self; mb_keys <+ "P_labile"; mb_values <+ P_labile;
			owners <+ self; mb_keys <+ "C_recalcitrant"; mb_values <+ C_recalcitrant;
			owners <+ self; mb_keys <+ "N_recalcitrant"; mb_values <+ N_recalcitrant;
			owners <+ self; mb_keys <+ "P_recalcitrant"; mb_values <+ P_recalcitrant;
		}
		ask Dam {
			owners <+ self; mb_keys <+ "C_dom (dom[2])"; mb_values <+ dom[2];
			owners <+ self; mb_keys <+ "N_dom (dom[0])"; mb_values <+ dom[0];
			owners <+ self; mb_keys <+ "P_dom (dom[1])"; mb_values <+ dom[1];
			owners <+ self; mb_keys <+ "N_dim (dim[0])"; mb_values <+ dim[0];
			owners <+ self; mb_keys <+ "P_dim (dim[1])"; mb_values <+ dim[1];
		}
		ask mb_microbes() {
			owners <+ self; mb_keys <+ "C (structural)"; mb_values <+ C;
			owners <+ self; mb_keys <+ "N (structural)"; mb_values <+ N;
			owners <+ self; mb_keys <+ "P (structural)"; mb_values <+ P;
			owners <+ self; mb_keys <+ "C_cytosol"; mb_values <+ cytosol_C;
			owners <+ self; mb_keys <+ "N_cytosol"; mb_values <+ cytosol_N;
			owners <+ self; mb_keys <+ "P_cytosol"; mb_values <+ cytosol_P;
		}
		ask Nematode {
			owners <+ self; mb_keys <+ "C_stomack"; mb_values <+ stomack_C;
			owners <+ self; mb_keys <+ "N_stomack"; mb_values <+ stomack_N;
			owners <+ self; mb_keys <+ "P_stomack"; mb_values <+ stomack_P;
		}
		ask MineralParticle {
			owners <+ self; mb_keys <+ "N"; mb_values <+ N;
			owners <+ self; mb_keys <+ "P"; mb_values <+ P;
		}

		mb_negative_count <- 0;
		// Garde-fou : "loop from: 0 to: -1" ne doit pas être exécuté si la liste est vide
		if (empty(mb_values)) { return; }
		loop i from: 0 to: length(mb_values) - 1 {
			string k <- mb_keys[i];
			float tol <- (k starts_with "C") ? tol_C : ((k starts_with "N") ? tol_N : tol_P);
			if (mb_values[i] < -tol) {
				mb_negative_count <- mb_negative_count + 1;
				if (mb_warning_count < mb_max_warnings) {
					write "[Positivité] Stock négatif au cycle " + cycle + " : " + owners[i].name + "." + k
						+ " = " + (mb_values[i] / #gram) + " g";
					mb_warning_count <- mb_warning_count + 1;
				}
			}
		}
		mb_negative_total <- mb_negative_total + mb_negative_count;
	}
}
