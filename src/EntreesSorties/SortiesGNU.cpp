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

#include "SortiesGNU.h"
#include "../Run.h"

using namespace std;
using namespace tinyxml2;

//***********************************************************************

SortiesGNU::SortiesGNU(){}

//***********************************************************************

SortiesGNU::SortiesGNU(string casTest, string run, XMLElement *element, string nomFichier, Entrees *entree) :
  Sorties(casTest, run, element, nomFichier, entree)
{
  m_nomFichierVisu = "visualisation.gnu";
}

//***********************************************************************

SortiesGNU::~SortiesGNU(){}

//***********************************************************************

void SortiesGNU::ecritSolutionSpecifique(Maillage *maillage, std::vector<Cellule *> *cellulesLvl)
{
  try {
    ofstream fluxFichier;
    string fichier = m_dossierSortie + creationNomFichierGNU(m_nomFichierResultats.c_str(), -1, rang, m_numFichier);
    fluxFichier.open(fichier.c_str());
    if (!fluxFichier) { throw ErreurECOGEN("Impossible d ouvrir le fichier " + fichier, __FILE__, __LINE__); }
    maillage->ecritSolutionGnuplot(cellulesLvl, fluxFichier);
    fluxFichier << endl;
    fluxFichier.close();

    //Creation du fichier gnuplot pour visualisation des resultats
    if (rang == 0) ecritScriptGnuplot(maillage->getGeometrie());
  }
  catch (ErreurECOGEN &) { throw; }
}

//*********************************************************************** 

void SortiesGNU::ecritScriptGnuplot(const int &dim)
{
  try {
    ofstream fluxFichier;

    fluxFichier.open((m_dossierSortie + m_nomFichierVisu).c_str());
    if (!fluxFichier) { throw ErreurECOGEN("Impossible d ouvrir le fichier " + m_dossierSortie + m_nomFichierVisu, __FILE__, __LINE__); }
    fluxFichier << "reset" << endl << "set style data lines" << endl << endl;

    int indice = 1 + dim; //premiere colonne fichier resultat == premiere donnee apres X

    //Gestion 2D/3D
    if (dim == 2) {
      fluxFichier << "set surface" << endl;
      fluxFichier << "set dgrid3d 50,50" << endl;
      fluxFichier << "set contour base" << endl;
      fluxFichier << "set cntrparam levels 25" << endl;
      fluxFichier << "show contour" << endl;
      //fluxFichier << "set view 180,180,1,1" << endl;

    }
    else if (dim == 3) { throw ErreurECOGEN("SortiesGNU::ecritScriptGnuplot : ecriture script gnuplot non prevu en 3D", __FILE__, __LINE__); }

    //1) Variables des phases
    //-----------------------
    for (int phase = 0; phase < m_run->getNombrePhases(); phase++)
    {
      //Variables scalaires
      for (int var = 1; var <= m_celluleRef.getPhase(phase)->getNombreScalaires(); var++) {
        fluxFichier << "set title '" << m_celluleRef.getPhase(phase)->renvoieNomScalaire(var) << "_" << m_celluleRef.getPhase(phase)->getEos()->retourneNom() << "'" << endl;
        ecritureBlocGnuplot(fluxFichier, indice, dim);
      } //Fin var scalaire
      //Variables vectorielles (u)
      for (int var = 1; var <= m_celluleRef.getPhase(phase)->getNombreVecteurs(); var++) {
        fluxFichier << "set title '" << m_celluleRef.getPhase(phase)->renvoieNomVecteur(var) << "_" << m_celluleRef.getPhase(phase)->getEos()->retourneNom() << "'" << endl;
        ecritureBlocGnuplot(fluxFichier, indice, dim);
      } //Fin var vectorielle
    } //Fin phase

   //2) Variables melange
   //--------------------
   //Variables scalaires
    for (int var = 1; var <= m_celluleRef.getMelange()->getNombreScalaires(); var++) {
      fluxFichier << "set title '" << m_celluleRef.getMelange()->renvoieNomScalaire(var) << "'" << endl;
      ecritureBlocGnuplot(fluxFichier, indice, dim);
    } //Fin var scalaire
    //Variables vectorielle
    for (int var = 1; var <= m_celluleRef.getMelange()->getNombreVecteurs(); var++) {
      fluxFichier << "set title '" << m_celluleRef.getMelange()->renvoieNomVecteur(var) << "'" << endl;
      ecritureBlocGnuplot(fluxFichier, indice, dim);
    } //Fin var vectorielle

    //3) Variables transports
    //-----------------------
    for (int var = 1; var <= m_celluleRef.getNombreTransports(); var++) {
      fluxFichier << "set title 'Transport" << var << "'" << endl;
      ecritureBlocGnuplot(fluxFichier, indice, dim);
    } //Fin var scalaire

    //4) Ecriture niveaux AMR
    //-----------------------
    fluxFichier << "set title 'Niveau AMR'" << endl;
    ecritureBlocGnuplot(fluxFichier, indice, dim);

    //5) Ecriture variable detection gradients
    //----------------------------------------
    fluxFichier << "set title 'Xi'" << endl;
    ecritureBlocGnuplot(fluxFichier, indice, dim);

    fluxFichier.close();

  }
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************

//***********************************************************************

void SortiesGNU::ecritureBlocGnuplot(std::ofstream &fluxFichier, int &indice, const int &dim)
{
  try {
    if (dim == 1) { fluxFichier << "plot"; }
    else { fluxFichier << "splot"; }
    for (int t = 0; t <= m_numFichier; t++) {
      for (int p = 0; p < Ncpu; p++) {
        fluxFichier << " \"" << creationNomFichierGNU(m_nomFichierResultats.c_str(), -1, p, t) << "\"";
        if (dim == 1) { fluxFichier << " u 1:" << indice; }
        else { fluxFichier << " u 1:2:" << indice; }
        if (p < Ncpu - 1 || t != m_numFichier) fluxFichier << ",\\";
        fluxFichier << endl;
      }
    }
    fluxFichier << "pause(-1)" << endl << endl;
    indice++;
  }
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************

string SortiesGNU::creationNomFichierGNU(const char* nom, int lvl, int proc, int numFichier, string nomVariable) const
{
  try {
    stringstream num;

    if (m_donneesSeparees) {
      throw ErreurECOGEN("SortiesGNU::creationNomFichierGNU : donnees Separees non prevu", __FILE__, __LINE__);
      return 0;
    }
    num << nom;
    //Gestion nomVariable
    if (nomVariable != "defaut") num << "_" << nomVariable << "_";
    //Gestion binaire
    if (m_ecritBinaire)
    {
      //num << "B64";
      throw ErreurECOGEN("SortiesGNU::creationNomFichierGNU : donnees binaires non prevu", __FILE__, __LINE__);
      return 0;
    }
    //Gestion cpu
    if (proc > -1) num << "_CPU" << proc;
    //Gestion niveau AMR
    if (lvl != -1)
    {
      //num << "_AMR" << lvl;
      throw ErreurECOGEN("SortiesGNU::creationNomFichierGNU : donnees AMR par niveau non prevu", __FILE__, __LINE__);
      return 0;
    }
    //Gestion numero de fichier resultat
    if (numFichier != -1) num << "_TIME" << numFichier;
    //Gestion extension
    num << ".out";
    return num.str();
  }
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************