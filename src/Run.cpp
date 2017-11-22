//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names 
//  are listed in the copyright file included with this source 
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published 
//  by the Free Software Foundation, either version 3 of the License, 
//  or (at your option) any later version.
//  
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).  
//  If not, see <http://www.gnu.org/licenses/>.

#include "Run.h"
#include "Outils.h"

using namespace std;
using namespace tinyxml2;

//***********************************************************************

Run::Run(std::string nomCasTest, const int &numero) : m_nombreTransports(0), m_evaporation(0), m_repriseFichier(0), m_dt(1e10), m_tempsPhys(0.), m_tempsCalc(0), m_iteration(0),
  m_casTest(nomCasTest), m_numTest(numero)
{}

//***********************************************************************

Run::~Run(){}

//***********************************************************************

void Run::initialisation(int argc, char* argv[])
{
  int eos_count(0); //Variable pour affecter un numero aux EOS

  //1) Initialisation parallele (a laisser meme en mono CPU)
  //--------------------------------------------------------
  Calcul_Parallele.initialisation(argc, argv);
  if (rang == 0){
    m_tempsInitial = clock();
  }
  if (Ncpu > 1){
    MPI_Barrier(MPI_COMM_WORLD);
    if (rang == 0) cout << "T" << m_numTest << " | Nombre de processeurs : " << Ncpu << endl;
  }

  //2) Lecture des fichiers d entrees au format XML
  //-----------------------------------------------
  vector<DomaineGeometrique*> domaines;
  vector<CondLim*> condLim;
  try {
    m_entree = new Entrees(this);
    m_entree->lectureEntreesXML(domaines, condLim);
  }
  catch (ErreurXML &) { throw; }
  BO = new Outils(m_nombrePhases);

  //3) Initialisation des donnees de maillage
  //-----------------------------------------
  m_maillage->attributLimites(condLim);
  try {
    m_dimension = m_maillage->initialiseGeometrie(&m_cellules, &m_bords, m_pretraitementParallele, m_ordre);
  }
  catch (ErreurECOGEN &) { throw; }
  int nombreCellules = m_maillage->getNombreCellules();
  int nombreCellulesTotales(nombreCellules);
  if (Ncpu>1) nombreCellulesTotales = m_maillage->getNombreCellulesTotal();

  //4) Initialisation des tableaux principaux selon modele et nombre de phases
  //--------------------------------------------------------------------------
  for (int i = 0; i < nombreCellulesTotales; i++) { m_cellules[i]->alloue(m_nombrePhases, m_nombreTransports, m_physAdd, m_modele); }
  //Attribution modele et pentes aux faces
  int nombreFaces = m_maillage->getNombreFaces();
  for (int i = 0; i < nombreFaces; i++) { m_bords[i]->associeModele(m_modele); }

  //5) Initialisation Physique : Remplissage donnees physiques
  //----------------------------------------------------------
  for (int i = 0; i < nombreCellulesTotales; i++) {
    m_cellules[i]->rempli(domaines);
  }
  try {
    if (m_repriseFichier) this->repriseFichier(m_iteration, m_dt, m_tempsPhys, m_tempsCalc); //Possibilite de reprise depuis fichier resultat
  }
  catch (ErreurECOGEN &) { throw; }
  // Remplissage du tableau eos
  m_cellules[0]->alloueEosYk(m_nombrePhases, m_modele);
  //Calcul ici les proprietes etendues des cellules (vitesses du son, energies, grandeurs de melange, etc.)
  for (int i = 0; i < nombreCellules; i++) {
    m_cellules[i]->calculsEtendus(m_nombrePhases);
  }

  //6) Allocation des Pentes et Cellules tampons pour problemes de Riemann
  //----------------------------------------------------------------------
  int allouePenteLocal = 0;
  for (int i = 0; i < nombreFaces; i++) { m_bords[i]->allouePentes(m_nombrePhases, m_nombreTransports, allouePenteLocal); }
  cellGauche = new Cellule; cellDroite = new Cellule;
  cellGauche->alloue(m_nombrePhases, m_nombreTransports, m_physAdd, m_modele);
  cellDroite->alloue(m_nombrePhases, m_nombreTransports, m_physAdd, m_modele);
  domaines[0]->rempli(cellGauche, m_nombrePhases, m_nombreTransports);
  domaines[0]->rempli(cellDroite, m_nombrePhases, m_nombreTransports);

  //7) Generation du(des) tableau(x) de cellules et bords non AMR ou AMR (pour chaque niveau)
  //-----------------------------------------------------------------------------------------
  m_maillage->genereTableauxCellulesBordsLvl(m_cellules, m_bords, &m_cellulesLvl, &m_bordsLvl);

  //8) Initialisation des communications persistantes pour parralele
  //----------------------------------------------------------------
	m_maillage->initialiseCommunicationsPersistantes(m_nombrePhases, m_nombreTransports, m_cellules, m_ordre);
  if (Ncpu > 1) {
    m_maillage->communicationsPrimitives(m_cellules, m_eos, 0);
    for (int i = 0; i < nombreCellules; i++) {
      m_cellules[i]->calculsEtendus(m_nombrePhases);
    }
  }
  
	//9) Initialisation Adaptive Mesh Refinement (AMR)
	//------------------------------------------------
	m_maillage->procedureRaffinementInitialisation(m_cellulesLvl, m_bordsLvl, m_physAdd, m_modele, m_nbMaillesTotalAMR, domaines, m_cellules, m_eos);

  for (unsigned int d = 0; d < domaines.size(); d++) { delete domaines[d]; }

  //10) Preparation des fichiers de sorties
  //---------------------------------------
  try {
    m_sortie->prepareSorties(*cellGauche);
    m_sortie->prepareSortiesInfos();
    for (unsigned int c = 0; c < m_coupes.size(); c++) m_coupes[c]->prepareSorties(*cellGauche);
	  if (rang == 0) m_sortie->ecritInfos();
	  m_sortie->ecritSolution(m_maillage, m_cellulesLvl);
	  for (unsigned int c = 0; c < m_coupes.size(); c++) m_coupes[c]->ecritSolution(m_maillage, m_cellulesLvl);
  }
  catch (ErreurXML &) { throw; }
  if (rang == 0) cout << " ...OK" << endl;

}

//***********************************************************************

void Run::repriseFichier(int &iteration, double &dt, double &tempsPhysique, clock_t &temps)
{
  //m_maillage->repriseCalcul(m_cellules, m_nombrePhases, m_nombreTransports, m_repriseFichier, m_ecritVTK, m_ecritXML, m_ecritBinaire);
  //FP//TODO// Reprise de calculs a faire
  
  int numFichier;
  ifstream fluxFichier;
  try {
    throw ErreurECOGEN("reprise impossible - non implementee'");
    fluxFichier.open("infosCalcul.out", ios::in);
    do {
      fluxFichier >> numFichier >> iteration >> tempsPhysique >> temps >> dt;
    } while (numFichier != m_repriseFichier && !fluxFichier.eof());
    if (fluxFichier.eof()) throw ErreurECOGEN("reprise impossible - verifier le fichier 'infosCalcul.out'");
  }
  catch (ErreurECOGEN &) { fluxFichier.close(); throw; }
  fluxFichier.close();
}

//***********************************************************************

void Run::resolution()
{
  int nombreFaces = m_maillage->getNombreFaces();
  int nbMaillesTotalAMRMax = m_nbMaillesTotalAMR;
  double dtMax;
  m_dt = 1.e-10;

  //---------------
  //Boucle en temps
  //---------------
  bool calculFini(false); bool ecriture(false);
  double ecritureSuivante(m_tempsPhys+m_freqTemps);
  while (!calculFini) {
    //Verification Erreurs
    try {
      this->verifieErreurs();
    }
    catch (ErreurECOGEN &) { throw; }
		
    //------------------- PROCEDURE INTEGRATION Run -------------------

    //Mise a zero des cons generale pour shema spatial sur dU/dt : permet de s affranchir du pas de temps
    for (unsigned int i = 0; i < m_cellulesLvl[0].size(); i++) { m_cellulesLvl[0][i]->miseAZeroConsGenerale(m_nombrePhases, m_nombreTransports); }
    //Procedure d integration Run
    dtMax = 1.e10;
    int lvlDep = 0;
    this->procedureIntegration(m_dt, lvlDep, dtMax, m_nbMaillesTotalAMR);

    //-------------------- CONTROLE ITERATIONS/TEMPS ---------------------

    //Still alive...
    if (m_iteration !=0 && m_iteration % 1000 == 0 && rang == 0) { cout << "Iteration " << m_iteration << " / Pas de temps " << m_dt << " / Avancement " << m_tempsPhys / m_tempsFinal*100. << "%" << endl; }

    m_tempsPhys += m_dt;
    m_iteration++;
    //Gestion Ecriture Resultats / Sortie boucle temporelle
    if (m_controleIterations) {
      if (m_iteration%m_freq == 0) { ecriture = true; }
      if (m_iteration >= m_nbIte) { calculFini = true; }
    }
    else {
      if (m_tempsPhys>ecritureSuivante) { ecriture = true; ecritureSuivante += m_freqTemps; }
      if (m_tempsPhys >= m_tempsFinal) { calculFini = true; }
    }

    //------------------------ ECRITURE FICHIERS -------------------------

    nbMaillesTotalAMRMax = max(nbMaillesTotalAMRMax, m_nbMaillesTotalAMR);
    if (ecriture) {
      clock_t tempsCalcul(clock() - m_tempsInitial);
      m_tempsInitial = clock();
      m_tempsCalc += tempsCalcul;
      //Ecritures generales
      if(rang==0) m_sortie->ecritInfos();
      m_sortie->ecritSolution(m_maillage, m_cellulesLvl);
      //Ecriture des coupes
      for (unsigned int c = 0; c < m_coupes.size(); c++) m_coupes[c]->ecritSolution(m_maillage, m_cellulesLvl);
      //Fin ecritures
      if (rang == 0) cout << " ...OK" << endl;
      ecriture = false;
    }

    //-------------------------- MISE A JOUR DT --------------------------

    //Mise a jour du pas de temps pour l iteration suivante
    m_dt = m_cfl * dtMax;
    if (Ncpu > 1) { Calcul_Parallele.calculDt(m_dt); }

  } //Fin boucle temps
  if (rang == 0) cout << "T" << m_numTest << " | ---------------------------------------" << endl;
  MPI_Barrier(MPI_COMM_WORLD);
  cout << "T" << m_numTest << " | Maximum de mailles sur CPU " << rang << " : " << nbMaillesTotalAMRMax << endl;
}

//***********************************************************************

void Run::procedureIntegration(double &dt, int lvl, double &dtMax, int &nbMaillesTotalAMR)
{
  //1) Calcul du pas de temps du niveau
  double dtLvl = dt * pow(2., -(double)lvl);

  //2) Procedure de (de)raffinement
  if (m_lvlMax > 0) { m_maillage->procedureRaffinement(m_cellulesLvl, m_bordsLvl, lvl, m_physAdd, m_modele, nbMaillesTotalAMR, m_cellules, m_eos); }

  //3) Calcul des pentes et gradients pour physiques additionnnelles
  //Fait ici pour avoir une mise a jour d'effectuer lors de l'execution de la procedure de niveau lvl+1 (donc pour les pentes plus besoin de les faire au debut de resolHyperboliqueO2)
  if (m_ordre == "ORDRE2") {
		for (unsigned int i = 0; i < m_bordsLvl[lvl].size(); i++) {
			if (!m_bordsLvl[lvl][i]->getSplit()) { m_bordsLvl[lvl][i]->calculPentes(m_nombrePhases, m_nombreTransports); }
		}
    if (Ncpu > 1) {
      m_maillage->communicationsPentes(m_cellules, lvl);
      if (lvl > 0) { m_maillage->communicationsPentes(m_cellules, lvl - 1); } //Comble un defaut d'un cas particulier non communique a temps (fantome niveau quelconque, cellule lvl l, cellule lvl l+1)
      //KS//FP// A reflechir si mieux a faire (a appliquer sur les 3 communications des pentes)
    }
  }
  if (lvl < m_lvlMax) {
    if (m_nombrePhysAdd) {
			for (unsigned int i = 0; i < m_cellulesLvl[lvl].size(); i++) {
				if (!m_cellulesLvl[lvl][i]->getSplit()) { m_cellulesLvl[lvl][i]->preparePhysAdd(); }
			}
    }

  //4) Appel a la procedure d'integration du niveau superieur
    this->procedureIntegration(dt, lvl + 1, dtMax, nbMaillesTotalAMR);
  }

  //5) Procedure d'avancement
  this->procedureAvancement(dtLvl, lvl, dtMax);

  //6) On refait une deuxieme fois l'enchainement pour les niveaux superieurs a 0
  if (lvl > 0) {
    if (m_ordre == "ORDRE2") {
			for (unsigned int i = 0; i < m_bordsLvl[lvl].size(); i++) {
        if (!m_bordsLvl[lvl][i]->getSplit()) { m_bordsLvl[lvl][i]->calculPentes(m_nombrePhases, m_nombreTransports); }
			}
      if (Ncpu > 1) {
        m_maillage->communicationsPentes(m_cellules, lvl);
        if (lvl > 0) { m_maillage->communicationsPentes(m_cellules, lvl - 1); }
      }
    }
    if (lvl < m_lvlMax) { this->procedureIntegration(dt, lvl + 1, dtMax, nbMaillesTotalAMR); }
    this->procedureAvancement(dtLvl, lvl, dtMax);
  }
}

//***********************************************************************

void Run::procedureAvancement(double &dt, int &lvl, double &dtMax) const
{
  //1) Resolution schema hyperbolique (Godunov ou MUSCL + Riemann)
  if (m_ordre == "ORDRE1") this->resolHyperbolique(dt, lvl, dtMax);
  else this->resolHyperboliqueO2(dt, lvl, dtMax);
  //2) Resolution volume fini des physiques ajoutees
  if (m_nombrePhysAdd) this->resolPhysiquesAdditionelles(dt, lvl);
  //3) Resolution/integration des termes sources sur modele non relaxe
  if (m_nombreSources) this->resolTermesSources(dt, lvl);
  //4) Resolution des relaxations eventuelles
  if (m_nombrePhases > 1) this->resolRelaxations(lvl); //FP//TODO// Gerer mieux les relaxations vis a vis des modeles
  //5) Moyenne des enfants dans la cellule mere
	if (lvl < m_lvlMax) { for (unsigned int i = 0; i < m_cellulesLvl[lvl].size(); i++) { m_cellulesLvl[lvl][i]->moyenneEnfantsDansParent(); } }
	//6) Communications finales
	if (Ncpu > 1) { m_maillage->communicationsPrimitives(m_cellules, m_eos, lvl); }
}

//***********************************************************************

void Run::resolHyperboliqueO2(double &dt, int &lvl, double &dtMax) const
{
  int nombreFaces = m_maillage->getNombreFaces();

  //1) Sauvegarde de m_cons pour combinaison Ordre2 et AMR
  //------------------------------------------------------
  for (unsigned int i = 0; i < m_cellulesLvl[lvl].size(); i++) { if (!m_cellulesLvl[lvl][i]->getSplit()) { m_cellulesLvl[lvl][i]->sauvegardeCons(m_nombrePhases, m_nombreTransports); } }

  //2) Schema semi-discret ordre2 spatial
  //-------------------------------------
  //Calcul de la somme des flux que l on stock dans m_cons de chaque cellule et determination du pas de temps hyperbolique max
  for (unsigned int i = 0; i < m_bordsLvl[lvl].size(); i++) { if (!m_bordsLvl[lvl][i]->getSplit()) { m_bordsLvl[lvl][i]->calculFlux(m_nombrePhases, m_nombreTransports, dtMax, *m_limiteurGlobal, *m_limiteurInterface); } }

  //3) Etape de prediction avec les pentes sur 1/2 pas de temps
  //-----------------------------------------------------------
  for (unsigned int i = 0; i < m_cellulesLvl[lvl].size(); i++) { if (!m_cellulesLvl[lvl][i]->getSplit()) { m_cellulesLvl[lvl][i]->predictionOrdre2(dt, m_nombrePhases, m_nombreTransports); } }

  //4) Recuperation de m_cons pour combinaison Ordre2 et AMR (remplace miseAZeroCons)
  //---------------------------------------------------------------------------------
  for (unsigned int i = 0; i < m_cellulesLvl[lvl].size(); i++) { if (!m_cellulesLvl[lvl][i]->getSplit()) { m_cellulesLvl[lvl][i]->recuperationCons(m_nombrePhases, m_nombreTransports); } }

  //5) Communication des vecPhasesO2
  //--------------------------------
  if (Ncpu > 1) {	m_maillage->communicationsPrimitives(m_cellules, m_eos, lvl, vecPhasesO2);	}

  //6) Nouveau calcul de pentes (optionnel mais stabilisant)
  //--------------------------------------------------------
  for (unsigned int i = 0; i < m_bordsLvl[lvl].size(); i++) { if (!m_bordsLvl[lvl][i]->getSplit()) { m_bordsLvl[lvl][i]->calculPentes(m_nombrePhases, m_nombreTransports, vecPhasesO2); } }
  if (Ncpu > 1) {
    m_maillage->communicationsPentes(m_cellules, lvl);
    if (lvl > 0) { m_maillage->communicationsPentes(m_cellules, lvl - 1); }
  }

  //7) Schema semi discret spatial sur les variables au demi pas de temps
  //---------------------------------------------------------------------
  //Calcul de la somme des flux que l on stock dans m_cons de chaque cellule et determination du pas de temps hyperbolique max
  for (unsigned int i = 0; i < m_bordsLvl[lvl].size(); i++) { if (!m_bordsLvl[lvl][i]->getSplit()) { m_bordsLvl[lvl][i]->calculFlux(m_nombrePhases, m_nombreTransports, dtMax, *m_limiteurGlobal, *m_limiteurInterface, vecPhasesO2); } }

  //8) Evolution temporelle
  //-----------------------
  for (unsigned int i = 0; i < m_cellulesLvl[lvl].size(); i++) {
    if (!m_cellulesLvl[lvl][i]->getSplit()) {
      m_cellulesLvl[lvl][i]->evolutionTemporelle(dt, m_nombrePhases, m_nombreTransports);   //Obtention des cons pour shema sur (Un+1-Un)/dt
      m_cellulesLvl[lvl][i]->construitPrim(m_nombrePhases);                                 //On peut reconstruire Prim a partir de m_cons
      m_cellulesLvl[lvl][i]->miseAZeroCons(m_nombrePhases, m_nombreTransports);             //Mise a zero des cons pour shema spatial sur dU/dt : permet de s affranchir du pas de temps
    }
  }
}

//***********************************************************************

void Run::resolHyperbolique(double &dt, int &lvl, double &dtMax) const
{
  //1) Schema semi-discret spatial
  //------------------------------
  //Calcul de la somme des flux que l on stock dans m_cons de chaque cellule et determination du pas de temps hyperbolique max
  for (unsigned int i = 0; i < m_bordsLvl[lvl].size(); i++) { if (!m_bordsLvl[lvl][i]->getSplit()) { m_bordsLvl[lvl][i]->calculFlux(m_nombrePhases, m_nombreTransports, dtMax, *m_limiteurGlobal, *m_limiteurInterface); } }

  //2) Evolution temporelle
  //-----------------------
  for (unsigned int i = 0; i < m_cellulesLvl[lvl].size(); i++) {
    if (!m_cellulesLvl[lvl][i]->getSplit()) {
      m_cellulesLvl[lvl][i]->evolutionTemporelle(dt, m_nombrePhases, m_nombreTransports);   //Obtention des cons pour shema sur (Un+1-Un)/dt
      m_cellulesLvl[lvl][i]->construitPrim(m_nombrePhases);                                 //On peut reconstruire Prim a partir de m_cons
      m_cellulesLvl[lvl][i]->miseAZeroCons(m_nombrePhases, m_nombreTransports);             //Mise a zero des cons pour shema spatial sur dU/dt : permet de s affranchir du pas de temps
    }
  }
}

//***********************************************************************

void Run::resolPhysiquesAdditionelles(double &dt, int &lvl) const
{
  //1) Preparation des grandeurs pour physiques additionelle (Calcul des gradients, etc) et communications
  //------------------------------------------------------------------------------------------------------
  if (Ncpu > 1) { m_maillage->communicationsPrimitives(m_cellules, m_eos, lvl); }
  for (unsigned int i = 0; i < m_cellulesLvl[lvl].size(); i++) { if (!m_cellulesLvl[lvl][i]->getSplit()) { m_cellulesLvl[lvl][i]->preparePhysAdd(); } }
  if (Ncpu > 1) { m_maillage->communicationsPhysAdd(m_physAdd, m_cellules, lvl); }

  //2) Calcul des flux des physiques additionnelles (Capillarite, Viscosite, Conductivite, ...)
  //-------------------------------------------------------------------------------------------
  //Calcul de la somme des flux des physiques additionnelles que l on stock dans m_cons de chaque cellule
  m_cellules[0]->miseAZeroFluxTemp(m_nombrePhases); // KS// Pas sur que ce soit encore utile, a verifier plus precisement !!!!! .... (avec cas test sans AMR)
  for (unsigned int pa = 0; pa < m_physAdd.size(); pa++) {
    for (unsigned int i = 0; i < m_bordsLvl[lvl].size(); i++) { if (!m_bordsLvl[lvl][i]->getSplit()) { m_bordsLvl[lvl][i]->calculFluxPhysAdd(m_nombrePhases, *m_physAdd[pa]); } }
    for (unsigned int i = 0; i < m_cellulesLvl[lvl].size(); i++) { if (!m_cellulesLvl[lvl][i]->getSplit()) { m_cellulesLvl[lvl][i]->ajoutNonConsPhysAdd(m_nombrePhases, *m_physAdd[pa]); } }
  }

  //3) Schema temporel physiques additionelles
  //------------------------------------------
  for (unsigned int i = 0; i < m_cellulesLvl[lvl].size(); i++) {
    if (!m_cellulesLvl[lvl][i]->getSplit()) {
      m_cellulesLvl[lvl][i]->evolutionTemporelle(dt, m_nombrePhases, m_nombreTransports);   //Obtention des cons pour shema sur (Un+1-Un)/dt
      m_cellulesLvl[lvl][i]->construitPrim(m_nombrePhases);                                 //On peut reconstruire Prim a partir de m_cons
      m_cellulesLvl[lvl][i]->miseAZeroCons(m_nombrePhases, m_nombreTransports);             //Mise a zero des cons pour shema spatial sur dU/dt : permet de s affranchir du pas de temps
    }
  }
}

//***********************************************************************

void Run::resolTermesSources(double &dt, int &lvl) const
{
  for (unsigned int i = 0; i < m_cellulesLvl[lvl].size(); i++) {
    if (!m_cellulesLvl[lvl][i]->getSplit()) {
      for (unsigned int s = 0; s < m_sources.size(); s++) { m_sources[s]->integreTermeSource(m_cellulesLvl[lvl][i], m_nombrePhases, dt); }
      m_cellulesLvl[lvl][i]->miseAZeroCons(m_nombrePhases, m_nombreTransports);
    }
  }
}

//***********************************************************************

void Run::resolRelaxations(int &lvl) const
{
  //Relaxation des pressions
  for (unsigned int i = 0; i < m_cellulesLvl[lvl].size(); i++) { if (!m_cellulesLvl[lvl][i]->getSplit()) { m_cellulesLvl[lvl][i]->relaxPressions(m_nombrePhases); } }
  //Re-initialisation de la fonction couleur (transport) avec alpha
  for (unsigned int pa = 0; pa < m_physAdd.size(); pa++) { m_physAdd[pa]->reinitialiseFonctionCouleur(m_cellulesLvl, lvl); }
  if (Ncpu > 1) { m_maillage->communicationsTransports(m_cellules, lvl); }
  for (unsigned int i = 0; i < m_cellulesLvl[lvl].size(); i++) { if (!m_cellulesLvl[lvl][i]->getSplit()) { m_cellulesLvl[lvl][i]->preparePhysAdd(); } }
  //Correction des energies et autres relaxations
  for (unsigned int i = 0; i < m_cellulesLvl[lvl].size(); i++) {
    if (!m_cellulesLvl[lvl][i]->getSplit()) {
      m_cellulesLvl[lvl][i]->correctionEnergie(m_nombrePhases);            //Correction des energies
      if (m_evaporation) m_cellulesLvl[lvl][i]->relaxPTMu(m_nombrePhases); //Relaxation des pressions, temperatures et potentiels chimiques
    }
  }
}

//***********************************************************************

void Run::verifieErreurs() const
{
  try {
    if (Ncpu > 1) {
      Calcul_Parallele.verifieEtatCPUs();
    }
    else if (erreurs.size() != 0) {
      throw ErreurECOGEN("Arret du code apres erreur... A voir");
    }
  }
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************

void Run::finalise()
{
  //Desallocations generales
  for (int i = 0; i < m_maillage->getNombreFaces(); i++) { delete m_bords[i]; } delete[] m_bords;
  for (int i = 0; i < m_maillage->getNombreCellulesTotal(); i++) { delete m_cellules[i]; } delete[] m_cellules;
  for (int i = 0; i < m_nombreEos; i++) { delete m_eos[i]; } delete[] m_eos;
  //Desallocations physiques additionnelles
  for (unsigned int pa = 0; pa < m_physAdd.size(); pa++) { delete m_physAdd[pa]; }
  for (unsigned int s = 0; s < m_sources.size(); s++) { delete m_sources[s]; }
  //Desallocations ordre 2
  if (m_ordre == "ORDRE2") {
    for (int k = 0; k < m_nombrePhases; k++) { delete pentesPhasesLocal1[k]; }
    for (int k = 0; k < m_nombrePhases; k++) { delete pentesPhasesLocal2[k]; }
    delete[] pentesPhasesLocal1;
    delete[] pentesPhasesLocal2;
    delete pentesMelangeLocal1;
    delete pentesMelangeLocal2;
    delete[] pentesTransportLocal1;
    delete[] pentesTransportLocal2;
  }
  //Desallocations parallele
	m_maillage->finaliseParallele(m_lvlMax);
  //Desallocations autres
  delete BO;
  delete cellGauche; delete cellDroite;
  delete m_maillage;
  delete m_modele;
  delete m_limiteurGlobal; delete m_limiteurInterface;
  delete m_entree; 
  delete m_sortie;
  delete[] GlMel_Yk;
  for (unsigned int s = 0; s < m_coupes.size(); s++) { delete m_coupes[s]; }
  //Desallocations AMR
  delete[] m_cellulesLvl;
  delete[] m_bordsLvl;
}

//***********************************************************************

void Run::arretApresErreur(void)
{
  cout << "EXCEPTION DETECTEE : La mort du programme a ete inconditionnelement brutale et sans aucune equivoque." << endl;
  exit(1);
}

//***********************************************************************

int Run::getNombrePhases() const { return m_nombrePhases; }

//***********************************************************************
