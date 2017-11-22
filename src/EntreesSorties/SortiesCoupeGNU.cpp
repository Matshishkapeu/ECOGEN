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

#include "SortiesCoupeGNU.h"

using namespace std;
using namespace tinyxml2;

//***************************************************************

SortiesCoupeGNU::SortiesCoupeGNU(){
  delete m_objet;
}

//***************************************************************

SortiesCoupeGNU::SortiesCoupeGNU(string casTest, string run, XMLElement *element, string nomFichier, TypeOG type, Entrees *entree)
{
  try {
    //Modification des attributs
    m_casTest = casTest;
    if (type == DROITE) { m_nomFichierResultats = "coupe1D"; }
    else if (type == PLAN) { m_nomFichierResultats = "coupe2D"; }
    else { throw ErreurECOGEN("SortiesCoupeGNU::SortiesCoupeGNU : type de coupe inconnu", __FILE__, __LINE__); }
    m_nomFichierVisu = "visualisation" + m_nomFichierResultats + ".gnu";
    m_dossierSortie = "./Resultats/" + run + "/coupes/";
    m_donneesSeparees = 0;
    m_numFichier = 0;
    m_entree = entree;
    m_run = m_entree->getRun();

    XMLElement *sousElement;
    XMLError erreur;

    double donnee;
    Coord point, vecteur;

    //Recuperation donnees de coupe 1D
    sousElement = element->FirstChildElement("point");
    if (sousElement == NULL) throw ErreurXMLElement("point", nomFichier, __FILE__, __LINE__);
    erreur = sousElement->QueryDoubleAttribute("x", &donnee); point.setX(donnee);
    if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("x", nomFichier, __FILE__, __LINE__);
    erreur = sousElement->QueryDoubleAttribute("y", &donnee); point.setY(donnee);
    if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("y", nomFichier, __FILE__, __LINE__);
    erreur = sousElement->QueryDoubleAttribute("z", &donnee); point.setZ(donnee);
    if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("z", nomFichier, __FILE__, __LINE__);
    if (type == DROITE) sousElement = element->FirstChildElement("vecDir");
    else if (type == PLAN) { sousElement = element->FirstChildElement("vecNormal"); }
    else { throw ErreurECOGEN("SortiesCoupeGNU::SortiesCoupeGNU : type de coupe inconnu", __FILE__, __LINE__); }
    if (sousElement == NULL) throw ErreurXMLElement("vecDir ou vecNormal", nomFichier, __FILE__, __LINE__);
    erreur = sousElement->QueryDoubleAttribute("x", &donnee); vecteur.setX(donnee);
    if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("x", nomFichier, __FILE__, __LINE__);
    erreur = sousElement->QueryDoubleAttribute("y", &donnee); vecteur.setY(donnee);
    if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("y", nomFichier, __FILE__, __LINE__);
    erreur = sousElement->QueryDoubleAttribute("z", &donnee); vecteur.setZ(donnee);
    if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("z", nomFichier, __FILE__, __LINE__);

    if (type == DROITE) { m_objet = new OGDroite(point, vecteur); }
    else if (type == PLAN) { m_objet = new OGPlan(point, vecteur); }
    else { throw ErreurECOGEN("SortiesCoupeGNU::SortiesCoupeGNU : type de coupe inconnu", __FILE__, __LINE__); }

  }
  catch (ErreurECOGEN &) { throw; }
}

//***************************************************************

SortiesCoupeGNU::~SortiesCoupeGNU(){}

//***********************************************************************

void SortiesCoupeGNU::ecritSolutionSpecifique(Maillage *maillage, std::vector<Cellule *> *cellulesLvl)
{
  ofstream fluxFichier;
  string fichier = m_dossierSortie + creationNomFichierGNU(m_nomFichierResultats.c_str(), -1, rang, m_numFichier);
  fluxFichier.open(fichier.c_str());
  maillage->ecritSolutionGnuplot(cellulesLvl, fluxFichier, m_objet);
  fluxFichier << endl;
  fluxFichier.close();

  try {
    //Creation du fichier gnuplot pour visualisation des resultats
    if (rang == 0) {
      if (m_objet->getType() == DROITE) { ecritScriptGnuplot(1); }
      else if (m_objet->getType() == PLAN) { ecritScriptGnuplot(2); }
      else { throw ErreurECOGEN("SortiesCoupeGNU::ecritSolutionSpecifique : type de coupe inconnu", __FILE__, __LINE__); }
    }
  }
  catch (ErreurECOGEN &) { throw; }
}

//***************************************************************
