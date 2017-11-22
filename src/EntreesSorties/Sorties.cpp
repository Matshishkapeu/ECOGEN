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

#include "Sorties.h"
#include "../Run.h"

using namespace std;
using namespace tinyxml2;

//***********************************************************************

Sorties::Sorties(){}

//***************************************************************
//Constructeur sortie a partir d une lecture au format XML modeSortie
//ex :	<modeSortie format="XML" binaire="false"/>

Sorties::Sorties(string casTest, string nomRun, XMLElement *element, string nomFichier, Entrees *entree) :
  m_casTest(casTest), m_dossierSortie(nomRun), m_numFichier(0), m_donneesSeparees(0), m_entree(entree)
{
  //Affectation pointeur run
  m_run = m_entree->getRun();

  //Noms communs
  //------------
  m_infosCalcul = "infoCalcul.out";
  m_infoMailles = "infosMaillage";
  m_nomFichierResultats = "result";
  m_nomFichierCollection = "collection";

  m_dossierSortie = "./Resultats/" + m_dossierSortie + "/";
  m_dossierSauvegardesEntrees = m_dossierSortie + "sauvegardesEntrees/";
  m_dossierSauvegardesInfosMailles = m_dossierSortie + "infosMaillage/";
  m_dossierCoupes = m_dossierSortie + "coupes/";

  //XMLElement *elementCoupe;
  XMLError erreur;

  //Recuperation mode Ecriture
  erreur = element->QueryBoolAttribute("binaire", &m_ecritBinaire);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("binaire", nomFichier, __FILE__, __LINE__);

  //Creation du dossier de sortie ou vidange /Macro selon OS Windows ou Linux
  if (rang == 0) {
    //Macro pour les interaction systeme (creation/destruction repertoires)
    #ifdef WIN32
      _mkdir("./Resultats");
      _mkdir(m_dossierSortie.c_str());
      _mkdir(m_dossierSauvegardesEntrees.c_str());
      _mkdir(m_dossierSauvegardesInfosMailles.c_str());
      _mkdir(m_dossierCoupes.c_str());
    #else
      mkdir("./Resultats", S_IRWXU);
      mkdir(m_dossierSortie.c_str(), S_IRWXU);
      mkdir(m_dossierSauvegardesEntrees.c_str(), S_IRWXU);
      mkdir(m_dossierSauvegardesInfosMailles.c_str(), S_IRWXU);
      mkdir(m_dossierCoupes.c_str(), S_IRWXU);
    #endif
    try {
      //Sauvegarde des fichiers d entrees
      IO::copieFichier(m_entree->getMain(), m_casTest, m_dossierSauvegardesEntrees);
      IO::copieFichier(m_entree->getMaillage(), m_casTest, m_dossierSauvegardesEntrees);
      IO::copieFichier(m_entree->getCI(), m_casTest, m_dossierSauvegardesEntrees);
      IO::copieFichier(m_entree->getModele(), m_casTest, m_dossierSauvegardesEntrees);
    }
    catch (ErreurECOGEN &) { throw; }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //Determination du mode Little / Big Endian
  //-----------------------------------------
  int entierTest = 42; //En binaire 0x2a
  char *chaineTest = reinterpret_cast<char*>(&entierTest);
  m_endianMode = "LittleEndian";
  if (chaineTest[0] != 0x2a) { m_endianMode = "BigEndian"; }
}

//***********************************************************************

Sorties::~Sorties(){}

//***********************************************************************

void Sorties::prepareSorties(const Cellule &cell)
{
  //Preparation de la cellule de reference
  //--------------------------------------
  m_celluleRef.alloue(m_run->m_nombrePhases, m_run->m_nombreTransports, m_run->m_physAdd, m_run->m_modele);
  for (int k = 0; k < m_run->m_nombrePhases; k++) { m_celluleRef.copiePhase(k, cell.getPhase(k)); }
  m_celluleRef.copieMelange(cell.getMelange());
  for (int k = 0; k < m_run->m_nombreTransports; k++) { m_celluleRef.setTransport(cell.getTransport(k).getValeur(), k); }

  //Preparation propres au type de sortie
  //-------------------------------------
  try {
    this->prepareSortieSpecifique();
  }
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************

void Sorties::prepareSortiesInfos()
{
  try {
    ofstream fluxFichier;
    //Fichier infosCalcul
    if (rang == 0) fluxFichier.open((m_dossierSortie + m_infosCalcul).c_str()); fluxFichier.close();
    //Fichiers infosMaillages
    string fichier = m_dossierSauvegardesInfosMailles + creationNomFichier(m_infoMailles.c_str(), -1, rang);
    fluxFichier.open(fichier.c_str()); fluxFichier.close();
  }
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************

void Sorties::ecritSolution(Maillage* maillage, vector<Cellule *> *cellulesLvl)
{
  try {
    //Ecriture specifiques a la sortie
    //--------------------------------
    this->ecritSolutionSpecifique(maillage, cellulesLvl);
  }
  catch (ErreurECOGEN &) { throw; }

  //Finalisation ecriture
  //---------------------
  MPI_Barrier(MPI_COMM_WORLD);
  m_numFichier++;
}

//***********************************************************************

void Sorties::ecritInfos()
{
  if (m_run->m_iteration > 0) {
    afficheInfoEcriture();
    sauvegardeInfos();
    sauvegardeInfosMailles();
  }
  cout << "T" << m_run->m_numTest << " | ecriture en cours fichier : " << m_numFichier << "... ";
}

//***********************************************************************

void Sorties::ecritJeuDonnees(std::vector<double> jeuDonnees, std::ofstream &fluxFichier, TypeDonnee typeDonnee)
{
  if (!m_ecritBinaire) {
    for (unsigned int k = 0; k < jeuDonnees.size(); k++) { fluxFichier << jeuDonnees[k] << " "; }
  }
  else {
    int donneeInt; float donneeFloat; double donneeDouble; char donneeChar;
    int taille;
    switch (typeDonnee) {
    case DOUBLE:
      taille = jeuDonnees.size()*sizeof(double); break;
    case FLOAT:
      taille = jeuDonnees.size()*sizeof(float); break;
    case INT:
      taille = jeuDonnees.size()*sizeof(int); break;
    case CHAR:
      taille = jeuDonnees.size()*sizeof(char); break;
    }
    IO::writeb64(fluxFichier, taille);
    char *chaineTampon = new char[taille]; int indice = 0;
    for (unsigned int k = 0; k < jeuDonnees.size(); k++) { 
      switch (typeDonnee) {
      case DOUBLE:
        donneeDouble = static_cast<double>(jeuDonnees[k]);
        IO::ajouteAlaChaine(chaineTampon, indice, donneeDouble);
        break;
      case FLOAT:
        donneeFloat = static_cast<float>(jeuDonnees[k]);
        IO::ajouteAlaChaine(chaineTampon, indice, donneeFloat);
        break;
      case INT:
        donneeInt = static_cast<int>(jeuDonnees[k]);
        IO::ajouteAlaChaine(chaineTampon, indice, donneeInt);
        break;
      case CHAR:
        donneeChar = static_cast<char>(jeuDonnees[k]);
        IO::ajouteAlaChaine(chaineTampon, indice, donneeChar);
        break;
      }
    }
    IO::writeb64Chaine(fluxFichier, chaineTampon, taille);
    delete[]chaineTampon;
  }
}

//***********************************************************************

void Sorties::afficheInfoEcriture() const
{
  double convDouble = static_cast<double>(m_run->m_tempsCalc) / CLOCKS_PER_SEC;
  int tempsCalcul = static_cast<int>(convDouble);

  cout << "T" << m_run->m_numTest << " | ---------------------------------------" << endl;
  cout << "T" << m_run->m_numTest << " | FICHIER RESULTAT : " << m_numFichier << ",  ITERATION " << m_run->m_iteration << endl;
  cout << "T" << m_run->m_numTest << " |     Temps physique    = " << m_run->m_tempsPhys << " s " << endl;
  cout << "T" << m_run->m_numTest << " |     Pas de temps      = " << m_run->m_dt << " s " << endl;
  //Affichage temps ecoule
  int seconde(tempsCalcul);
  if (seconde < 60)
  {
    cout << "T" << m_run->m_numTest << " |     Temps ecoule     = " << seconde << " s " << endl;
  }
  else
  {
    int minute(seconde / 60);
    seconde = seconde % 60;
    if (minute <60)
    {
      cout << "T" << m_run->m_numTest << " |     Temps ecoule     = " << minute << " min " << seconde << " s " << endl;
    }
    else
    {
      int heure(minute / 60);
      minute = minute % 60;
      cout << "T" << m_run->m_numTest << " |     Temps ecoule     = " << heure << " h " << minute << " min " << seconde << " s " << endl;
    }
  }
  //Estimation temps restant
  //A faire...
  cout << "T" << m_run->m_numTest << " | ---------------------------------------" << endl;
}

//***********************************************************************

void Sorties::sauvegardeInfos() const
{
  ofstream fluxFichier;
  if (rang == 0) {
    fluxFichier.open((m_dossierSortie + m_infosCalcul).c_str(), ios::app);
    fluxFichier << m_numFichier << " " << m_run->m_iteration << " " << m_run->m_tempsPhys << " " << m_run->m_tempsCalc << " " << m_run->m_dt << endl;
    fluxFichier.close();
  }
}

//***********************************************************************

void Sorties::sauvegardeInfosMailles() const
{
  try {
    ofstream fluxFichier;
    string fichier = m_dossierSauvegardesInfosMailles + creationNomFichier(m_infoMailles.c_str(), -1, rang);
    fluxFichier.open(fichier.c_str(), ios::app);
    fluxFichier << m_run->m_nbMaillesTotalAMR << endl;
    fluxFichier.close();
  }
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************

string Sorties::creationNomFichier(const char* nom, int lvl, int proc, int numFichier) const
{

  try {
    stringstream num;
    num << nom;
    //Gestion cpu
    if (proc > -1) num << "_CPU" << proc;
    //Gestion niveau AMR
    if (lvl != -1)
    { throw ErreurECOGEN("Sorties::creationNomFichier : donnees AMR par niveau non prevu", __FILE__, __LINE__); }
    //Gestion numero de fichier resultat
    if (numFichier != -1) num << "_TIME" << numFichier;
    //Gestion extension
    num << ".out";
    return num.str();
  }
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************