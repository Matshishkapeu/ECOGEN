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

#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "MaillageNonStruct.h"
#include "../Erreurs.h"

using namespace std;
using namespace tinyxml2;

//***********************************************************************

MaillageNonStruct::MaillageNonStruct(const string &fichierMaillage) :
  Maillage(),
  m_fichierMaillage(fichierMaillage),
  m_nomMaillage(fichierMaillage),
  m_nombreNoeuds(0),
  m_nombreNoeudsInternes(0),
  m_nombreElementsInternes(0),
  m_nombreElementsFantomes(0),
  m_nombreElementsCommunicants(0),
  m_nombreCellulesFantomes(0),
  m_nombreFacesParallele(0),
  m_nombreElements0D(0),
  m_nombreElements1D(0),
  m_nombreElements2D(0),
  m_nombreElements3D(0),
  m_nombreSegments(0),
  m_nombreTriangles(0),
  m_nombreQuadrangles(0),
  m_nombreTetraedres(0),
  m_nombrePyramides(0),
  m_nombrePoints(0),
  m_nombreHexaedres(0)
{
  m_nomMaillage = m_fichierMaillage;
  m_nomMaillage.resize(m_nomMaillage.size() - 4); //On enleve l extension
  m_type = UNS;
}

//***********************************************************************

MaillageNonStruct::~MaillageNonStruct(){
  for (int f = 0; f < m_nombreFacesTotal; f++) { delete m_faces[f]; }
  delete[] m_faces;
  for (int e = 0; e < m_nombreElements; e++) { delete m_elements[e]; }
  delete[] m_elements;
  delete[] m_noeuds;
  for (unsigned int l = 0; l < m_lim.size(); l++) { delete m_lim[l]; }
}

//***********************************************************************

void MaillageNonStruct::attributLimites(std::vector<CondLim*> &condLim)
{
  //Recherche numero physique de limite le plus grand
  int maxNumLim(0);
  for (unsigned int i = 0; i < condLim.size(); i++) {
    if (condLim[i]->getNumPhys() > maxNumLim) maxNumLim = condLim[i]->getNumPhys();
  }
  //Attribution des limites dans le tableau m_lim dans l'ordre
  int limite(1);
  int limiteTrouvee(0);
  while (limite <= maxNumLim) {
    for (unsigned int i = 0; i < condLim.size(); i++) {
      if (condLim[i]->getNumPhys() == limite) {
        m_lim.push_back(condLim[i]);
        limiteTrouvee = 1; break;
      }
    }
    if (!limiteTrouvee) m_lim.push_back(new CondLimAbs(limite));
    limite++; limiteTrouvee = 0;
  }
}

//***********************************************************************

int MaillageNonStruct::initialiseGeometrie(Cellule ***cellules, BordDeMaille ***bord, bool pretraitementParallele, string ordreCalcul)
{
  try {
    if (Ncpu == 1) { this->initialiseGeometrieMonoCPU(cellules, bord, ordreCalcul); }
    else {
      //Pretraitement du fichier de maillage par le CPU 0
      if (pretraitementParallele) {
        if (rang == 0) { this->pretraitementFichierMaillageGmsh(); }
        MPI_Barrier(MPI_COMM_WORLD);
      }
      this->initialiseGeometrieParallele(cellules, bord, ordreCalcul);
    }
    return m_geometrie;
  }
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************

void MaillageNonStruct::initialiseGeometrieMonoCPU(Cellule ***cellules, BordDeMaille ***bord, string ordreCalcul)
{
  try {
    //1) Lecture noeuds et elements
    //-----------------------------
    vector<ElementNS*>* voisinsNoeuds; //Dimensions : nb_noeuds
    this->lectureGeometrieGmsh(&voisinsNoeuds);  //Remplissage m_noeuds et m_elements

    cout << "------------------------------------------------------" << endl;
    cout << " B) CONSTRUCTION DE LA GEOMETRIE EN COURS..." << endl;

    //2) Attribution des cellules a leur element geometrique
    //------------------------------------------------------
    //Comptage mailles et limites
    if (m_nombreElements3D == 0 && m_nombreElements2D == 0) //Cas 1D
    {
      m_nombreCellulesCalcul = m_nombreElements1D;
      m_nombreFacesLimites = m_nombreElements0D;
      m_geometrie = 1;
    }
    else if (m_nombreElements3D == 0) //Cas 2D
    {
      m_nombreCellulesCalcul = m_nombreElements2D;
      m_nombreFacesLimites = m_nombreElements1D;
      m_geometrie = 2;
    }
    else //Cas 3D
    {
      m_nombreCellulesCalcul = m_nombreElements3D;
      m_nombreFacesLimites = m_nombreElements2D;
      m_geometrie = 3;
    }
    m_nombreCellulesTotal = m_nombreCellulesCalcul;

    //Dimensionnement du tableau de cellules
    (*cellules) = new Cellule*[m_nombreCellulesCalcul];

    //Attribution Absorption aux limites manquantes
    unsigned int nbLimites(0);
    for (int i = 0; i < m_nombreFacesLimites; i++) {
      unsigned int appPhys(m_elements[i]->getAppartenancePhysique());
      if (appPhys > nbLimites) { nbLimites = appPhys; }
    }
    for (unsigned int i = m_lim.size(); i < nbLimites; i++) {
      m_lim.push_back(new CondLimAbs);
    }

    //Attribution des cellules aux elements et comptage faces internes
    m_nombreFacesInternes = 0;
    for (int i = 0; i < m_nombreCellulesCalcul; i++)
    {
      if (ordreCalcul == "ORDRE1") { (*cellules)[i] = new Cellule; }
      else { (*cellules)[i] = new CelluleO2; }
      (*cellules)[i]->setElement(m_elements[i + m_nombreFacesLimites], i);
      m_nombreFacesInternes += m_elements[i + m_nombreFacesLimites]->getNombreFaces();
    }

    //-----------Ajout construction voisins des mailles pour ordre 2----------
    vector<ElementNS *>* voisins; //vecteur temporaire des voisins
    voisins = new vector<ElementNS *>[m_nombreCellulesCalcul];
    for (int i = 0; i < m_nombreCellulesCalcul; i++)
    {
      ElementNS *e(m_elements[i + m_nombreFacesLimites]); //choix element
      //cout << "el " << i + m_nombreFacesLimites << " : " ;
      //vector<ElementNS *> voisins; //vecteur temporaire des voisins
      //1) Construction du vecteur de voisins
      for (int n = 0; n < e->getNombreNoeuds(); n++) { //boucle noeud de element e
        int numNoeud = e->getNumNoeud(n);
        for (unsigned int v = 0; v < voisinsNoeuds[numNoeud].size(); v++) { //boucle voisin du noeud n
          bool ajoute(true);
          if (voisinsNoeuds[numNoeud][v]->getIndice() == e->getIndice()) ajoute = false;
          //else if(voisinsNoeuds[numNoeud][v]->getIndice() < m_nombreFacesLimites) ajoute = false;
          else {
            for (unsigned int vo = 0; vo < voisins[i].size(); vo++) {
              if (voisinsNoeuds[numNoeud][v]->getIndice() == voisins[i][vo]->getIndice()) {
                ajoute = false; break;
              }
            }
          }
          if (ajoute) {
            voisins[i].push_back(voisinsNoeuds[numNoeud][v]);
            //cout << voisinsNoeuds[numNoeud][v]->getIndice() << " ";
          }
        }
      }
      //cout << endl;
    }
    //------------------------------------------------------------------------

    m_nombreFacesInternes -= m_nombreFacesLimites; //On enleve les limites
    m_nombreFacesInternes /= 2; //Les faces internes sont toute comptees 2 fois => on retabli la verite !
    m_nombreFacesTotal = m_nombreFacesInternes + m_nombreFacesLimites;

    double volTot(0.);
    for (int i = 0; i < m_nombreCellulesCalcul; i++)
    {
      volTot += (*cellules)[i]->getElement()->getVolume();
    }
    cout << "Volume total : " << volTot << endl;

    //3) Construction de la table de connectivite
    //-------------------------------------------
    //Dimensionnement du tableau de faces
    (*bord) = new BordDeMaille*[m_nombreFacesTotal];
    m_faces = new FaceNS*[m_nombreFacesTotal];
    int **facesTemp; int *sommeNoeudsTemp; //On cre un tableau temporaire de faces pour accelerer la recherche d'existance
    facesTemp = new int*[m_nombreFacesTotal + 1];
    sommeNoeudsTemp = new int[m_nombreFacesTotal + 1];
    //Determination du nombre de noeuds max pour les faces
    int tailleFace; //Sera initialise a la taille maximale
    if (m_nombreElements3D != 0)
    {
      if (m_nombreQuadrangles != 0) { tailleFace = 4; }
      else if (m_nombreTriangles != 0) { tailleFace = 3; }
      else { Erreurs::messageErreur("Probleme dans initialiseGeometrieMonoCPU pour initialisation du tableau facesTemp"); }
    }
    else if (m_nombreElements2D != 0) { tailleFace = 2; }
    else { tailleFace = 1; }
    for (int i = 0; i < m_nombreFacesTotal + 1; i++)
    {
      facesTemp[i] = new int[tailleFace];
    }

    //Faces internes
    int indiceMaxFaces(0);
    clock_t tTemp(clock());
    float t1(0.);
    cout << "  1/Construction des faces ..." << endl;
    int frequenceImpression(max((m_nombreElements - m_nombreFacesLimites) / 10, 1));
    for (int i = m_nombreFacesLimites; i < m_nombreElements; i++)
    {
      if ((i - m_nombreFacesLimites) % frequenceImpression == 0) { cout << "    " << (100 * (i - m_nombreFacesLimites) / (m_nombreElements - m_nombreFacesLimites)) << "% ... " << endl; }
      m_elements[i]->construitFaces(m_noeuds, m_faces, indiceMaxFaces, facesTemp, sommeNoeudsTemp);
    }
    for (int i = 0; i < m_nombreFacesTotal + 1; i++) { delete facesTemp[i]; }
    delete[] facesTemp; delete[] sommeNoeudsTemp;
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    cout << "    OK en " << t1 << " secondes" << endl;

    //Limites
    cout << "  2/Attribution des elements limites aux faces limites ..." << endl;
    tTemp = clock();
    frequenceImpression = max(m_nombreFacesLimites / 10, 1);
    for (int i = 0; i < m_nombreFacesLimites; i++)
    {
      if (i%frequenceImpression == 0)
      {
        cout << "    " << (100 * i / m_nombreFacesLimites) << "% ... " << endl;// << " -> "; 
      }
      //Attribution de la limite
      m_elements[i]->attributFaceLimite(m_noeuds, m_faces, indiceMaxFaces);
    }
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    cout << "    OK en " << t1 << " secondes" << endl;

    //Liaison Geometrie/Bords de calcul
    cout << "  3/Liaisons Geometriques -> Physiques ..." << endl;
    tTemp = clock();
    int iMailleG, iMailleD;
    for (int i = 0; i < m_nombreFacesTotal; i++)
    {
      if (m_faces[i]->getEstLimite())
      {
        int appPhys(m_faces[i]->getElementDroite()->getAppartenancePhysique() - 1); //appartenance - 1 pour tableau commencant a zero
        if (appPhys >= static_cast<int>(m_lim.size()) || appPhys < 0) { Erreurs::messageErreur("Nombre de conditions aux limites non adapte"); }
        m_lim[appPhys]->creeLimite(&(*bord)[i]);
        (*bord)[i]->setFace(m_faces[i]);
        iMailleG = m_faces[i]->getElementGauche()->getIndice() - m_nombreFacesLimites;
        iMailleD = iMailleG;
      }
      else
      {
        if (ordreCalcul == "ORDRE1") { (*bord)[i] = new BordDeMaille; }
        else { (*bord)[i] = new BordDeMailleO2; }
        (*bord)[i]->setFace(m_faces[i]);
        iMailleG = m_faces[i]->getElementGauche()->getIndice() - m_nombreFacesLimites;
        iMailleD = m_faces[i]->getElementDroite()->getIndice() - m_nombreFacesLimites;
      }
      (*bord)[i]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
      (*cellules)[iMailleG]->ajouteBord((*bord)[i]);
      (*cellules)[iMailleD]->ajouteBord((*bord)[i]);

      //-------------Ici il faut trouver les mailles voisines pour le calcul ordre 2 multipentes ------------
      if (ordreCalcul != "ORDRE1") {

        //GAUCHE) Recherche des mailles a gauche
        //**************************************
        //ElementNS *eG(m_faces[i]->getElementGauche()); //Element a gauche
        //int indice(eG->getIndice() - m_nombreFacesLimites);


        //MaillageNonStruct::rechercheElementsArrieres(eG,m_faces[i], (*bord)[i], voisins[indice], *cellules);
        //MaillageNonStruct::rechercheElementsAvants(eG, m_faces[i], (*bord)[i], voisins[indice], *cellules);
        //
        //
        ////A enlever a terme...
        //double cos(0.),a,b,beta1,beta2;
        //ElementNS *eT;
        //Coord vG = eG->vecteur(m_faces[i]);
        //Coord MB1, MB2, MB;



        //Coord vG = eG->vecteur(m_faces[i]);
        //ElementNS *eT; //Element a tester
        ////1) Recherche du premier element arriere gauche BG1M
        ////---------------------------------------------------
        //double cos(0.);
        //for (unsigned int v = 0; v < voisins[indice].size(); v++) {
        //  eT = voisins[indice][v];
        //  Coord vT = eT->vecteur(eG);
        //  double cosTemp(Coord::cos(vT, vG));
        //  if (cosTemp >= cos) {
        //    cos = cosTemp;
        //    if (eT->getIndice() < m_nombreFacesLimites) { //Si c est une limite on remet BG1M a NULL
        //      (*bord)[i]->setB(BG1M, 0);
        //    }
        //    else { //Sinon on met a jour avec la maille la plus loin
        //      Cellule *c((*cellules)[eT->getNumCelluleAssociee()]);
        //      (*bord)[i]->setB(BG1M, c);
        //    }
        //  }
        //}
        ////2) Recherche du second element arriere gauche BG2M
        ////---------------------------------------------------
        //cos = 0.;
        //if ((*bord)[i]->getB(BG1M) != 0) {  //Cas non bord
        //  Element *e1((*bord)[i]->getB(BG1M)->getElement());
        //  Coord v1 = e1->vecteur(eG);
        //  Coord sin1(Coord::sin(v1, vG));
        //  if (fabs(sin1.norme()) <= 1e-8) {
        //    (*bord)[i]->setB(BG2M, (*bord)[i]->getB(BG1M)); // Si le cos est nul, meme cellule pour B1 et B2
        //  }
        //  else {
        //    for (unsigned int v = 0; v < voisins[indice].size(); v++) {
        //      eT = voisins[indice][v];
        //      if (eT == e1) continue; // voisin suivant si voisin deja utilise 
        //      Coord vT = eT->vecteur(eG);
        //      Coord sinTemp(Coord::sin(vT, vG));

        //      if (sinTemp.scalaire(sin1) <= 0.) {
        //        double cosTemp(Coord::cos(vT, vG));
        //        if (cosTemp >= cos) {
        //          cos = cosTemp;
        //          if (eT->getIndice() < m_nombreFacesLimites) { //Si c est une limite on remet BG1M a NULL
        //            (*bord)[i]->setB(BG2M, 0);
        //          }
        //          else {  //Sinon on met a jour avec la 2 eme maille la plus loin
        //            Cellule *c((*cellules)[eT->getNumCelluleAssociee()]);
        //            (*bord)[i]->setB(BG2M, c);
        //          }
        //        }
        //      } //fin sin*sin <0
        //    }  //fin boucle voisins
        //  } //fin if 
        //} //Fin if bord

        ////Determination des ponderations arrieres
        //Coord MB1; if ((*bord)[i]->getB(BG1M) != 0) { MB1.creeVecteur(m_faces[i]->getPos(), (*bord)[i]->getB(BG1M)->getPosition()); }
        //Coord MB2; if ((*bord)[i]->getB(BG2M) != 0) { MB2.creeVecteur(m_faces[i]->getPos(), (*bord)[i]->getB(BG2M)->getPosition()); }
        //Coord MB; MB.creeVecteur(m_faces[i]->getPos(), eG->getPosition());
        //double a, b, beta1(1.), beta2(0.);
        //a = (MB1.vectoriel(MB)).norme() / MB.norme();
        //b = (MB2.vectoriel(MB)).norme() / MB.norme();
        //if (fabs(a + b) > 1e-8) { beta1 = b / (a + b); beta2 = 1. - beta1; }
        //(*bord)[i]->setBeta(betaG1M, beta1);
        //(*bord)[i]->setBeta(betaG2M, beta2);
        ////Calcul de la distance du point M (bord de maille) au point H (barycentre)
        //if ((*bord)[i]->getB(BG1M) != 0) {
        //  cos = Coord::cos(MB, MB1);
        //  double d, e, c;
        //  d = cos*MB1.norme();
        //  Coord B1B2; B1B2 = MB1 - MB2;
        //  e = B1B2.norme()*beta1;
        //  c = sqrt(e*e - a*a);
        //  d += c;
        //  (*bord)[i]->setDistanceH(distanceHGM, d);
        //}
        ////cout << eG->getIndice() << " : " << endl;
        ////if ((*bord)[i]->getB(BG1M) != 0) cout << "BG1M = " << (*bord)[i]->getB(BG1M)->getElement()->getIndice() << " " << (*bord)[i]->getDistanceH(distanceHGM) << endl;
        ////if ((*bord)[i]->getB(BG2M) != 0) cout << "BG2M = " << (*bord)[i]->getB(BG2M)->getElement()->getIndice() << " " << (*bord)[i]->getDistanceH(distanceHGM) << endl;

        ////3) Recherche du premier element avant gauche BG1P
        ////-------------------------------------------------
        //cos = 0.;
        //for (unsigned int v = 0; v < voisins[indice].size(); v++) {
        //  eT = voisins[indice][v];
        //  Coord vT = eT->vecteur(eG);
        //  double cosTemp(Coord::cos(vT, vG));
        //  if (cosTemp <= cos) {
        //    cos = cosTemp;
        //    if (eT->getIndice() < m_nombreFacesLimites) { //Si c est une limite on remet BG1P a NULL
        //      (*bord)[i]->setB(BG1P, 0);
        //    }
        //    else { //Sinon on met a jour avec la maille la plus loin
        //      Cellule *c((*cellules)[eT->getNumCelluleAssociee()]);
        //      (*bord)[i]->setB(BG1P, c);
        //    }
        //  }
        //}
        ////4) Recherche du second element avant gauche BG2P
        ////------------------------------------------------
        //cos = 0.;
        //if ((*bord)[i]->getB(BG1P) != 0) {  //Cas non bord
        //  Element *e1((*bord)[i]->getB(BG1P)->getElement());
        //  Coord v1 = e1->vecteur(eG);
        //  Coord sin1(Coord::sin(v1, vG));
        //  if (fabs(sin1.norme()) <= 1e-8) {
        //    (*bord)[i]->setB(BG2P, (*bord)[i]->getB(BG1P)); // Si le cos est nul, meme cellule pour B1 et B2
        //  }
        //  else {
        //    for (unsigned int v = 0; v < voisins[indice].size(); v++) {
        //      eT = voisins[indice][v];
        //      if (eT == e1) continue; // voisin suivant si voisin deja utilise 
        //      Coord vT = eT->vecteur(eG);
        //      Coord sinTemp(Coord::sin(vT, vG));

        //      if (sinTemp.scalaire(sin1) <= 0.) {
        //        double cosTemp(Coord::cos(vT, vG));
        //        if (cosTemp <= cos) {
        //          cos = cosTemp;
        //          if (eT->getIndice() < m_nombreFacesLimites) { //Si c est une limite on remet BG2P a NULL
        //            (*bord)[i]->setB(BG2P, 0);
        //          }
        //          else {  //Sinon on met a jour avec la 2 eme maille la plus loin
        //            Cellule *c((*cellules)[eT->getNumCelluleAssociee()]);
        //            (*bord)[i]->setB(BG2P, c);
        //          }
        //        }
        //      } //fin sin*sin <0
        //    }  //fin boucle voisins
        //  } //fin if 
        //} //Fin if bord

        ////Determination des ponderations avant
        //MB1 = 0.; if ((*bord)[i]->getB(BG1P) != 0) { MB1.creeVecteur((*bord)[i]->getB(BG1P)->getPosition(), m_faces[i]->getPos()); }
        //MB2 = 0.; if ((*bord)[i]->getB(BG2P) != 0) { MB2.creeVecteur((*bord)[i]->getB(BG2P)->getPosition(), m_faces[i]->getPos()); }
        //MB.creeVecteur(eG->getPosition(), m_faces[i]->getPos());
        //beta1 = 1., beta2 = 0.;
        //a = (MB1.vectoriel(MB)).norme() / MB.norme();
        //b = (MB2.vectoriel(MB)).norme() / MB.norme();
        //if (fabs(a + b) > 1e-8) { beta1 = b / (a + b); beta2 = 1. - beta1; }
        //(*bord)[i]->setBeta(betaG1P, beta1);
        //(*bord)[i]->setBeta(betaG2P, beta2);
        ////Calcul de la distance du point M (bord de maille) au point H (barycentre)
        //if ((*bord)[i]->getB(BG1P) != 0) {
        //  cos = Coord::cos(MB, -1.*MB1);
        //  double d, e, c;
        //  d = cos*MB1.norme();
        //  Coord B1B2; B1B2 = MB1 - MB2;
        //  e = B1B2.norme()*beta1;
        //  c = sqrt(e*e - a*a);
        //  d += c;
        //  (*bord)[i]->setDistanceH(distanceHGP, d);
        //}
        ////cout << eG->getIndice() << " : " << endl;
        ////if ((*bord)[i]->getB(BG1P) != 0) cout << "BG1P = " << (*bord)[i]->getB(BG1P)->getElement()->getIndice() << " " << (*bord)[i]->getDistanceH(distanceHGP) << endl;
        ////if ((*bord)[i]->getB(BG2P) != 0) cout << "BG2P = " << (*bord)[i]->getB(BG2P)->getElement()->getIndice() << " " << (*bord)[i]->getDistanceH(distanceHGP) << endl;


        ////DROITE) Recherche des mailles a droite
        ////**************************************
        //ElementNS *eD(m_faces[i]->getElementDroite()); //Element a droite
        //indice = eD->getIndice() - m_nombreFacesLimites;
        //if (indice >= 0) {
        //  Coord vD = eD->vecteur(m_faces[i]);
        //  //1) Recherche du premier element arriere droite BD1M
        //  //---------------------------------------------------
        //  cos = 0.;
        //  for (unsigned int v = 0; v < voisins[indice].size(); v++) {
        //    eT = voisins[indice][v];
        //    Coord vT = eT->vecteur(eD);
        //    double cosTemp(Coord::cos(vT, vD));
        //    if (cosTemp >= cos) {
        //      cos = cosTemp;
        //      if (eT->getIndice() < m_nombreFacesLimites) { //Si c est une limite on remet BG1M a NULL
        //        (*bord)[i]->setB(BD1M, 0);
        //      }
        //      else { //Sinon on met a jour avec la maille la plus loin
        //        Cellule *c((*cellules)[eT->getNumCelluleAssociee()]);
        //        (*bord)[i]->setB(BD1M, c);
        //      }
        //    }
        //  }
        //  //2) Recherche du second element arriere droite BD2M
        //  //---------------------------------------------------
        //  cos = 0.;
        //  if ((*bord)[i]->getB(BD1M) != 0) {  //Cas non bord
        //    Element *e1((*bord)[i]->getB(BD1M)->getElement());
        //    Coord v1 = e1->vecteur(eD);
        //    Coord sin1(Coord::sin(v1, vD));
        //    if (fabs(sin1.norme()) <= 1e-8) {
        //      (*bord)[i]->setB(BD2M, (*bord)[i]->getB(BD1M)); // Si le cos est nul, meme cellule pour B1 et B2
        //    }
        //    else {
        //      for (unsigned int v = 0; v < voisins[indice].size(); v++) {
        //        eT = voisins[indice][v];
        //        if (eT == e1) continue; // voisin suivant si voisin deja utilise 
        //        Coord vT = eT->vecteur(eD);
        //        Coord sinTemp(Coord::sin(vT, vD));

        //        if (sinTemp.scalaire(sin1) <= 0.) {
        //          double cosTemp(Coord::cos(vT, vD));
        //          if (cosTemp >= cos) {
        //            cos = cosTemp;
        //            if (eT->getIndice() < m_nombreFacesLimites) { //Si c est une limite on remet BG1M a NULL
        //              (*bord)[i]->setB(BD2M, 0);
        //            }
        //            else {  //Sinon on met a jour avec la 2 eme maille la plus loin
        //              Cellule *c((*cellules)[eT->getNumCelluleAssociee()]);
        //              (*bord)[i]->setB(BD2M, c);
        //            }
        //          }
        //        } //fin sin*sin <0
        //      }  //fin boucle voisins
        //    } //fin if 
        //  } //Fin if bord

        //  //Determination des ponderations arrieres
        //  MB1 = 0.; if ((*bord)[i]->getB(BD1M) != 0) { MB1.creeVecteur(m_faces[i]->getPos(), (*bord)[i]->getB(BD1M)->getPosition()); }
        //  MB2 = 0.; if ((*bord)[i]->getB(BD2M) != 0) { MB2.creeVecteur(m_faces[i]->getPos(), (*bord)[i]->getB(BD2M)->getPosition()); }
        //  MB.creeVecteur(m_faces[i]->getPos(), eD->getPosition());
        //  beta1 = 1.; beta2 = 0.;
        //  a = (MB1.vectoriel(MB)).norme() / MB.norme();
        //  b = (MB2.vectoriel(MB)).norme() / MB.norme();
        //  if (fabs(a + b) > 1e-8) { beta1 = b / (a + b); beta2 = 1. - beta1; }
        //  (*bord)[i]->setBeta(betaD1M, beta1);
        //  (*bord)[i]->setBeta(betaD2M, beta2);
        //  //Calcul de la distance du point M (bord de maille) au point H (barycentre)
        //  if ((*bord)[i]->getB(BD1M) != 0) {
        //    cos = Coord::cos(MB, MB1);
        //    double d, e, c;
        //    d = cos*MB1.norme();
        //    Coord B1B2; B1B2 = MB1 - MB2;
        //    e = B1B2.norme()*beta1;
        //    c = sqrt(e*e - a*a);
        //    d += c;
        //    (*bord)[i]->setDistanceH(distanceHDM, d);
        //  }
        //  //cout << eD->getIndice() << " : " << endl;
        //  //if ((*bord)[i]->getB(BD1M) != 0) cout << "BD1M = " << (*bord)[i]->getB(BD1M)->getElement()->getIndice() << " " << (*bord)[i]->getDistanceH(distanceHDM) << endl;
        //  //if ((*bord)[i]->getB(BD2M) != 0) cout << "BD2M = " << (*bord)[i]->getB(BD2M)->getElement()->getIndice() << " " << (*bord)[i]->getDistanceH(distanceHDM) << endl;

        //  //3) Recherche du premier element avant droite BD1P
        //  //-------------------------------------------------
        //  cos = 0.;
        //  for (unsigned int v = 0; v < voisins[indice].size(); v++) {
        //    eT = voisins[indice][v];
        //    Coord vT = eT->vecteur(eD);
        //    double cosTemp(Coord::cos(vT, vD));
        //    if (cosTemp <= cos) {
        //      cos = cosTemp;
        //      if (eT->getIndice() < m_nombreFacesLimites) { //Si c est une limite on remet BG1M a NULL
        //        (*bord)[i]->setB(BD1P, 0);
        //      }
        //      else { //Sinon on met a jour avec la maille la plus loin
        //        Cellule *c((*cellules)[eT->getNumCelluleAssociee()]);
        //        (*bord)[i]->setB(BD1P, c);
        //      }
        //    }
        //  }
        //  //4) Recherche du second element avant droite BD2P
        //  //------------------------------------------------
        //  cos = 0.;
        //  if ((*bord)[i]->getB(BD1P) != 0) {  //Cas non bord
        //    Element *e1((*bord)[i]->getB(BD1P)->getElement());
        //    Coord v1 = e1->vecteur(eD);
        //    Coord sin1(Coord::sin(v1, vD));
        //    if (fabs(sin1.norme()) <= 1e-8) {
        //      (*bord)[i]->setB(BD2P, (*bord)[i]->getB(BD1P)); // Si le cos est nul, meme cellule pour B1 et B2
        //    }
        //    else {
        //      for (unsigned int v = 0; v < voisins[indice].size(); v++) {
        //        eT = voisins[indice][v];
        //        if (eT == e1) continue; // voisin suivant si voisin deja utilise 
        //        Coord vT = eT->vecteur(eD);
        //        Coord sinTemp(Coord::sin(vT, vD));

        //        if (sinTemp.scalaire(sin1) <= 0.) {
        //          double cosTemp(Coord::cos(vT, vD));
        //          if (cosTemp <= cos) {
        //            cos = cosTemp;
        //            if (eT->getIndice() < m_nombreFacesLimites) { //Si c est une limite on remet BG1M a NULL
        //              (*bord)[i]->setB(BD2P, 0);
        //            }
        //            else {  //Sinon on met a jour avec la 2 eme maille la plus loin
        //              Cellule *c((*cellules)[eT->getNumCelluleAssociee()]);
        //              (*bord)[i]->setB(BD2P, c);
        //            }
        //          }
        //        } //fin sin*sin <0
        //      }  //fin boucle voisins
        //    } //fin if 
        //  } //Fin if bord
        //} //Fin si limite

        ////Determination des ponderations avant
        //MB1 = 0.; if ((*bord)[i]->getB(BD1P) != 0) { MB1.creeVecteur((*bord)[i]->getB(BD1P)->getPosition(), m_faces[i]->getPos()); }
        //MB2 = 0.; if ((*bord)[i]->getB(BD2P) != 0) { MB2.creeVecteur((*bord)[i]->getB(BD2P)->getPosition(), m_faces[i]->getPos()); }
        //MB.creeVecteur(eD->getPosition(), m_faces[i]->getPos());
        //beta1 = 1., beta2 = 0.;
        //a = (MB1.vectoriel(MB)).norme() / MB.norme();
        //b = (MB2.vectoriel(MB)).norme() / MB.norme();
        //if (fabs(a + b) > 1e-8) { beta1 = b / (a + b); beta2 = 1. - beta1; }
        //(*bord)[i]->setBeta(betaD1P, beta1);
        //(*bord)[i]->setBeta(betaD2P, beta2);
        ////Calcul de la distance du point M (bord de maille) au point H (barycentre)
        //if ((*bord)[i]->getB(BD1P) != 0) {
        //  cos = Coord::cos(MB, -1.*MB1);
        //  double d, e, c;
        //  d = cos*MB1.norme();
        //  Coord B1B2; B1B2 = MB1 - MB2;
        //  e = B1B2.norme()*beta1;
        //  c = sqrt(e*e - a*a);
        //  d += c;
        //  (*bord)[i]->setDistanceH(distanceHDP, d);
        //}
        //cout << eD->getIndice() << " : " << endl;
        //if ((*bord)[i]->getB(BD1P) != 0) cout << "BD1P = " << (*bord)[i]->getB(BD1P)->getElement()->getIndice() << " " << (*bord)[i]->getDistanceH(distanceHDP) << endl;
        //if ((*bord)[i]->getB(BD2P) != 0) cout << "BD2P = " << (*bord)[i]->getB(BD2P)->getElement()->getIndice() << " " << (*bord)[i]->getDistanceH(distanceHDP) << endl;



        //cout << "******" << endl;
        //cout << eG->getIndice() << " " << eD->getIndice() << " : " << endl;
        //if ((*bord)[i]->getB(BG1M) != 0) cout << "BG1M = " << (*bord)[i]->getB(BG1M)->getElement()->getIndice() << endl;
        //if ((*bord)[i]->getB(BG2M) != 0) cout << "BG2M = " << (*bord)[i]->getB(BG2M)->getElement()->getIndice() << endl;
        //if ((*bord)[i]->getB(BG1P) != 0) cout << "BG1P = " << (*bord)[i]->getB(BG1P)->getElement()->getIndice() << endl;
        //if ((*bord)[i]->getB(BG2P) != 0) cout << "BG2P = " << (*bord)[i]->getB(BG2P)->getElement()->getIndice() << endl;
        //if ((*bord)[i]->getB(BD1M) != 0) cout << "BD1M = " << (*bord)[i]->getB(BD1M)->getElement()->getIndice() << endl;
        //if ((*bord)[i]->getB(BD2M) != 0) cout << "BD2M = " << (*bord)[i]->getB(BD2M)->getElement()->getIndice() << endl;
        //if ((*bord)[i]->getB(BD1P) != 0) cout << "BD1P = " << (*bord)[i]->getB(BD1P)->getElement()->getIndice() << endl;
        //if ((*bord)[i]->getB(BD2P) != 0) cout << "BD2P = " << (*bord)[i]->getB(BD2P)->getElement()->getIndice() << endl;

      } //Fin preparation voisins ordre 2
      //---------------------------------------------------------------------------------------

    } //Fin face
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    cout << "    OK en " << t1 << " secondes" << endl;
    cout << "... FIN CONSTRUCTION DE LA GEOMETRIE " << endl;
    cout << "------------------------------------------------------" << endl;

    delete[] voisins;
    delete[] voisinsNoeuds;
  }
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************

void MaillageNonStruct::initialiseGeometrieParallele(Cellule ***cellules, BordDeMaille ***bord, string ordreCalcul)
{

  clock_t tempsTotal(clock());

  //1) Lecture noeuds et elements
  //-----------------------------
  try {
    this->lectureGeometrieGmshParallele(); //Remplissage de m_noeuds et m_elements
    if (rang == 0)
    {
      cout << "------------------------------------------------------" << endl;
      cout << " B) CONSTRUCTION DE LA GEOMETRIE EN COURS..." << endl;
    }

    //2) Attribution des cellules a leur element geometrique
    //------------------------------------------------------
    //Comptage mailles et limites
    if (m_nombreElements3D == 0 && m_nombreElements2D == 0) //Cas 1D
    {
      m_nombreCellulesCalcul = m_nombreElements1D;
      m_nombreFacesLimites = m_nombreElements0D;
      m_geometrie = 1;
    }
    else if (m_nombreElements3D == 0) //Cas 2D
    {
      m_nombreCellulesCalcul = m_nombreElements2D;
      m_nombreFacesLimites = m_nombreElements1D;
      m_geometrie = 2;
    }
    else //Cas 3D
    {
      m_nombreCellulesCalcul = m_nombreElements3D;
      m_nombreFacesLimites = m_nombreElements2D;
      m_geometrie = 3;
    }

    //Dimensionnement du tableau de cellules
    m_nombreCellulesTotal = m_nombreCellulesCalcul + m_nombreElementsFantomes;
    (*cellules) = new Cellule*[m_nombreCellulesTotal];

    //Attribution Absorption aux limites manquantes
    unsigned int nbLimites(0);
    for (int i = 0; i < m_nombreFacesLimites; i++) {
      unsigned int appPhys(m_elements[i]->getAppartenancePhysique());
      if (appPhys > nbLimites) { nbLimites = appPhys; }
    }
    for (unsigned int i = m_lim.size(); i < nbLimites; i++) {
      m_lim.push_back(new CondLimAbs);
    }

    //Attribution des cellules aux elements et comptage faces internes
    m_nombreFacesInternes = 0;
    for (int i = 0; i < m_nombreCellulesTotal; i++)
    {
      if (ordreCalcul == "ORDRE1") { (*cellules)[i] = new Cellule; }
      else { (*cellules)[i] = new CelluleO2; }
      (*cellules)[i]->setElement(m_elements[i + m_nombreFacesLimites], i);
      if (i < m_nombreCellulesCalcul) { m_nombreFacesInternes += m_elements[i + m_nombreFacesLimites]->getNombreFaces(); }
    }
    m_nombreFacesInternes -= m_nombreFacesLimites + m_nombreFacesParallele; //On enleve les limites et les faces communicantes
    m_nombreFacesInternes /= 2; //Les faces internes sont toutes comptees 2 fois => on retabli la verite !
    m_nombreFacesTotal = m_nombreFacesInternes + m_nombreFacesLimites + m_nombreFacesParallele;

    //3) Construction de la table de connectivite interne
    //---------------------------------------------------
    //Dimensionnement des tableaux de faces
    (*bord) = new BordDeMaille*[m_nombreFacesTotal];
    m_faces = new FaceNS*[m_nombreFacesTotal];

    //On cree un tableau temporaire de faces pour accelerer la recherche d'existance
    int **facesTemp; int *sommeNoeudsTemp;
    facesTemp = new int*[m_nombreFacesTotal + 1];
    sommeNoeudsTemp = new int[m_nombreFacesTotal + 1];
    //Determination du nombre de noeuds max pour les faces
    int tailleFace; //Sera initialise a la taille maximale
    if (m_nombreElements3D != 0)
    {
      if (m_nombreQuadrangles != 0) { tailleFace = 4; }
      else if (m_nombreTriangles != 0) { tailleFace = 3; }
      else { Erreurs::messageErreur("Probleme dans initialiseGeometrieMonoCPU pour initialisation du tableau facesTemp"); }
    }
    else if (m_nombreElements2D != 0) { tailleFace = 2; }
    else { tailleFace = 1; }
    for (int i = 0; i < m_nombreFacesTotal + 1; i++) // Le +1 est utilise pour la recherche existance de face
    {
      facesTemp[i] = new int[tailleFace]; // Inconnue sur le nombre de points d une face (maxi 4 a priori)
    }

    //Faces internes
    //--------------
    int indiceMaxFaces(0);
    MPI_Barrier(MPI_COMM_WORLD);
    int frequenceImpression;
    clock_t tTemp(clock());
    float t1(0.);
    if (rang == 0)
    {
      cout << "  1/Construction des faces ..." << endl;
      frequenceImpression = max((m_nombreElementsInternes - m_nombreFacesLimites) / 10, 1);
    }
    for (int i = m_nombreFacesLimites; i < m_nombreElementsInternes; i++)
    {
      if (rang == 0 && (i - m_nombreFacesLimites) % frequenceImpression == 0) { cout << "    " << (100 * (i - m_nombreFacesLimites) / (m_nombreElementsInternes - m_nombreFacesLimites)) << "% ... " << endl; }
      //Construction
      m_elements[i]->construitFaces(m_noeuds, m_faces, indiceMaxFaces, facesTemp, sommeNoeudsTemp);
    }
    for (int i = 0; i < m_nombreFacesTotal + 1; i++) { delete facesTemp[i]; }
    delete[] facesTemp; delete[] sommeNoeudsTemp;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rang == 0)
    {
      tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
      cout << "    OK en " << t1 << " secondes" << endl;
    }

    //Limites
    //-------
    tTemp = clock();
    if (rang == 0)
    {
      cout << "  2/Attribution des elements limites aux faces limites ..." << endl;
      frequenceImpression = max(m_nombreFacesLimites / 10, 1);
    }
    for (int i = 0; i < m_nombreFacesLimites; i++)
    {
      if (rang == 0 && i%frequenceImpression == 0) { cout << "    " << (100 * i / m_nombreFacesLimites) << "% ... " << endl; }
      //Attribution de la limite
      m_elements[i]->attributFaceLimite(m_noeuds, m_faces, indiceMaxFaces);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rang == 0)
    {
      tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
      cout << "    OK en " << t1 << " secondes" << endl;
    }

    //Communications
    //--------------
    tTemp = clock();
    if (rang == 0)
    {
      cout << "  3/Attribution des ghost cells aux faces communicantes ..." << endl;
      frequenceImpression = max((m_nombreElements - m_nombreElementsInternes) / 10, 1);
    }
    for (int i = m_nombreElementsInternes; i < m_nombreElements; i++)
    {
      if (rang == 0 && (i - m_nombreElementsInternes) % frequenceImpression == 0) { cout << "    " << (100 * (i - m_nombreElementsInternes) / (m_nombreElements - m_nombreElementsInternes)) << "% ... " << endl; }
      //Attribution de la limite communicante
      m_elements[i]->attributFaceCommunicante(m_noeuds, m_faces, indiceMaxFaces, m_nombreNoeudsInternes);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rang == 0)
    {
      tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
      cout << "    OK en " << t1 << " secondes" << endl;
    }

    //Liaison Geometrie/Bords de calcul
    //---------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    tTemp = clock();
    if (rang == 0)
    {
      cout << "  4/Liaisons Geometriques -> Physiques ..." << endl;
      frequenceImpression = max(m_nombreFacesTotal / 10, 1);
    }
    int iMailleG, iMailleD;

    for (int i = 0; i < m_nombreFacesTotal; i++)
    {
      if (rang == 0 && i%frequenceImpression == 0) { cout << "    " << 100 * i / m_nombreFacesTotal << "% ... " << endl; }
      //Faces limites (limite physique ou communication)
      if (m_faces[i]->getEstLimite())
      {
        //Communication
        if (m_faces[i]->getEstComm())
        {
          if (ordreCalcul == "ORDRE1") { (*bord)[i] = new BordDeMaille; }
          else { (*bord)[i] = new BordDeMailleO2; }
          (*bord)[i]->setFace(m_faces[i]);
          iMailleG = m_faces[i]->getElementGauche()->getNumCelluleAssociee();
          iMailleD = m_faces[i]->getElementDroite()->getNumCelluleAssociee();
        }
        //Limite physique
        else
        {
          int appPhys(m_faces[i]->getElementDroite()->getAppartenancePhysique() - 1); //appartenance - 1 pour tableau commencant a zero
          if (appPhys >= static_cast<int>(m_lim.size()) || appPhys < 0) { throw ErreurECOGEN("Nombre de conditions aux limites non adapte", __FILE__, __LINE__); }
          m_lim[appPhys]->creeLimite(&(*bord)[i]);
          (*bord)[i]->setFace(m_faces[i]);
          iMailleG = m_faces[i]->getElementGauche()->getNumCelluleAssociee();
          iMailleD = iMailleG;
        }
      }
      //face interne au domaine
      else
      {
        if (ordreCalcul == "ORDRE1") { (*bord)[i] = new BordDeMaille; }
        else { (*bord)[i] = new BordDeMailleO2; }
        (*bord)[i]->setFace(m_faces[i]);
        iMailleG = m_faces[i]->getElementGauche()->getNumCelluleAssociee();
        iMailleD = m_faces[i]->getElementDroite()->getNumCelluleAssociee();
      }
      (*bord)[i]->initialise((*cellules)[iMailleG], (*cellules)[iMailleD]);
      (*cellules)[iMailleG]->ajouteBord((*bord)[i]);
      (*cellules)[iMailleD]->ajouteBord((*bord)[i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rang == 0)
    {
      tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
      cout << "    OK en " << t1 << " secondes" << endl;
    }

    //4) Construction de la table de connectivite parallele CPUs
    //----------------------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    tTemp = clock();
    if (rang == 0)
    {
      cout << "  5/Construction tables connectivite des CPUs ..." << endl;
      frequenceImpression = max((m_nombreElements - m_nombreFacesLimites) / 10, 1);
    }

    vector<int> nombreElementsAEnvoyer(Ncpu);
    for (int i = 0; i < Ncpu; i++) { nombreElementsAEnvoyer[i] = 0; }
    vector< vector<int> > elementsAEnvoyer(Ncpu);
    vector<int> nombreElementsARecevoir(Ncpu);
    for (int i = 0; i < Ncpu; i++) { nombreElementsARecevoir[i] = 0; }
    vector< vector<int> > elementsARecevoir(Ncpu);

    for (int i = m_nombreFacesLimites; i < m_nombreElements; i++)
    {
      int nombreAutresCPUs = m_elements[i]->getNombreAutresCPU();
      if (nombreAutresCPUs != 0) //Elements qui communique
      {
        if (m_elements[i]->getCPU() == rang) //L element appartient au CPU
        {
          for (int p = 0; p < nombreAutresCPUs; p++)
          {
            int numCPU = m_elements[i]->getAutreCPU(p);
            nombreElementsAEnvoyer[numCPU]++;
            elementsAEnvoyer[numCPU].push_back(i);
          }
        }
        else // L'element appartient a un autre CPU
        {
          int numCPU = m_elements[i]->getCPU();
          nombreElementsARecevoir[numCPU]++;
          elementsARecevoir[numCPU].push_back(i);
        }
      }
    }
    int *tampon;
    for (int v = 0; v < Ncpu; v++)
    {
      if (nombreElementsAEnvoyer[v] != 0) Calcul_Parallele.setVoisin(v);
      tampon = new int[nombreElementsAEnvoyer[v]];
      for (int i = 0; i < nombreElementsAEnvoyer[v]; i++) { tampon[i] = elementsAEnvoyer[v][i] - m_nombreFacesLimites; }
      Calcul_Parallele.setElementsAEnvoyer(v, tampon, nombreElementsAEnvoyer[v]);
      delete[] tampon;
      tampon = new int[nombreElementsARecevoir[v]];
      for (int i = 0; i < nombreElementsARecevoir[v]; i++) { tampon[i] = elementsARecevoir[v][i] - m_nombreFacesLimites; }
      Calcul_Parallele.setElementsARecevoir(v, tampon, nombreElementsARecevoir[v]);
      delete[] tampon;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rang == 0)
    {
      tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
      cout << "    OK en " << t1 << " secondes" << endl;
    }

    if (rang == 0)
    {
      cout << "... FIN CONSTRUCTION DE LA GEOMETRIE ";
      tempsTotal = clock() - tempsTotal; t1 = static_cast<float>(tempsTotal) / CLOCKS_PER_SEC;
      cout << " Temps construction geometrique total : " << t1 << " secondes" << endl;
      cout << "------------------------------------------------------" << endl;
    }
  }
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************

void MaillageNonStruct::pretraitementFichierMaillageGmsh()
{
  ElementNS **elementsGlobal;
  Coord *noeudsGlobal;
  int frequenceImpression(0);
  int nombreNoeudsGlobal(0), nombreElementsGlobal(0);
  int nombreElements0D(0), nombreElements1D(0), nombreElements2D(0), nombreElements3D(0);

  try {
    //Ouverture du fichier de maillage
    //-------------------------------------
    cout << "------------------------------------------------------" << endl;
    cout << " 0) PRETRAITEMENT DU FICHIER DE MAILLAGE " + m_fichierMaillage + " EN COURS..." << endl;
    clock_t tempsTotal(clock());
    m_fichierMaillage = "./libMaillages/" + m_fichierMaillage;
    ifstream fichierMaillage(m_fichierMaillage.c_str(), ios::in);
    if (!fichierMaillage) { throw ErreurECOGEN("fichier maillage absent :" + m_fichierMaillage, __FILE__, __LINE__); }
    string ligneCourante;
    getline(fichierMaillage, ligneCourante);
    getline(fichierMaillage, ligneCourante);
    getline(fichierMaillage, ligneCourante);
    getline(fichierMaillage, ligneCourante);

    //1) Stockage de la grille de point dans tableau temporaire de Noeuds
    //-------------------------------------------------------
    cout << "  1/lecture des noeuds du maillage ...";
    fichierMaillage >> nombreNoeudsGlobal;
    fichierMaillage.ignore(1, '\n');
    noeudsGlobal = new Coord[nombreNoeudsGlobal];
    int inutile(0); double x, y, z;
    for (int i = 0; i < nombreNoeudsGlobal; i++)
    {
      fichierMaillage >> inutile >> x >> y >> z;
      noeudsGlobal[i].setXYZ(x, y, z);
    }
    fichierMaillage.ignore(1, '\n');
    getline(fichierMaillage, ligneCourante);
    getline(fichierMaillage, ligneCourante);
    cout << "OK" << endl;

    //2) Recuperation des elements 1D/2D/3D dans le tableau temporaire d'elements et comptage
    //----------------------------------------------------------------------------
    cout << "  2/lecture des elements 1D/2D/3D ...";
    fichierMaillage >> nombreElementsGlobal;
    fichierMaillage.ignore(1, '\n');
    //Allocation tableau d elements
    elementsGlobal = new ElementNS*[nombreElementsGlobal];
    //Lecture des elements et attributions proprietes geometriques
    int posDebutElement = static_cast<int>(fichierMaillage.tellg()); //reperage debut des elements pour retour rapide
    for (int i = 0; i < nombreElementsGlobal; i++)
    {
      this->lectureElementGmsh(noeudsGlobal, fichierMaillage, &elementsGlobal[i]);
      //Comptage des elements
      if (elementsGlobal[i]->getTypeGmsh() == 15) { nombreElements0D++; }
      else if (elementsGlobal[i]->getTypeGmsh() == 1) { nombreElements1D++; }
      else if (elementsGlobal[i]->getTypeGmsh() <= 3) { nombreElements2D++; }
      else if (elementsGlobal[i]->getTypeGmsh() <= 7) { nombreElements3D++; }
      else { Erreurs::messageErreur("Type element du .msh non gere dans ECOGEN"); }
    }
    fichierMaillage.close();

    //Comptage mailles et limites
    int nombreElementsFrontiere;
    //Comptage mailles et limites
    if (nombreElements3D == 0 && nombreElements2D == 0) //Cas 1D
    {
      nombreElementsFrontiere = nombreElements0D;
    }
    else if (nombreElements3D == 0) //Cas 2D
    {
      nombreElementsFrontiere = nombreElements1D;
    }
    else //Cas 3D
    {
      nombreElementsFrontiere = nombreElements2D;
    }

    cout << "OK" << endl;
    cout << "  -----------------------------------" << endl;
    cout << "    INFORMATIONS GLOBALES MAILLAGE :" << endl;
    cout << "  -----------------------------------" << endl;
    cout << "    nombre de noeuds de maillage : " << nombreNoeudsGlobal << endl;
    cout << "    nombre elements  : " << nombreElementsGlobal << endl;

    //3)Reperage et attribution par CPU
    //---------------------------------
    cout << "  3/Attribution des noeuds par CPU ..." << endl;
    clock_t tTemp(clock()); float t1(0);
    frequenceImpression = max(nombreElementsGlobal / 10, 1);
    vector< vector<int> > noeudsCPU(Ncpu);
    vector< vector<int> > elementsCPU(Ncpu);
    int numCPU; int noeudCourant; bool noeudExiste;
    int numCPUMax(0); //pour verification maillage adapte ou non
    //Recherche des noeuds formant le maillage interne (hors fantomes) pour chaque CPU
    for (int i = 0; i < nombreElementsGlobal; i++)
    {
      if (i%frequenceImpression == 0) { cout << "    " << 100 * i / nombreElementsGlobal << "% ... " << endl; }
      numCPU = elementsGlobal[i]->getCPU();
      if (numCPU > numCPUMax) numCPUMax = numCPU;
      elementsCPU[numCPU].push_back(i); //Remplissage numero des elements propres au CPU
      for (int n = 0; n < elementsGlobal[i]->getNombreNoeuds(); n++)
      {
        noeudCourant = elementsGlobal[i]->getNumNoeud(n);
        noeudExiste = false;
        for (unsigned int j = 0; j < noeudsCPU[numCPU].size(); j++)
        {
          if (noeudsCPU[numCPU][j] == noeudCourant) { noeudExiste = true; break; }
        }
        if (!noeudExiste) { noeudsCPU[numCPU].push_back(noeudCourant); }
      }
    }
    int *nombreNoeudsInterne = new int[Ncpu];
    for (int p = 0; p < Ncpu; p++)
    {
      nombreNoeudsInterne[p] = noeudsCPU[p].size();
    }
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    cout << "    OK en " << t1 << " secondes" << endl;

    //Verification adaptation au maillage
    if (numCPUMax != Ncpu - 1) throw ErreurECOGEN("fichier maillage .msh non adapte au nombre de CPU - Generer le maillage et relancer le test", __FILE__, __LINE__);

    //4) Creation Tableau de faces
    //-----------------------------
    //Test estimation nombre de faces (approximative)
    cout << "  4/Creation tableau de faces ..." << endl;
    tTemp = clock();
    frequenceImpression = max((nombreElementsGlobal - nombreElementsFrontiere) / 10, 1);
    int *nombreFacesTemp = new int[Ncpu];
    int *iMaxFaces = new int[Ncpu];
    for (int p = 0; p < Ncpu; p++)
    {
      iMaxFaces[p] = 0;
      nombreFacesTemp[p] = 0;
      for (unsigned int i = 0; i < elementsCPU[numCPU].size(); i++)
      {
        nombreFacesTemp[p] += elementsGlobal[elementsCPU[numCPU][i]]->getNombreFaces();
      }
    }
    int ***facesTemp2 = new int**[Ncpu];
    int **sommeNoeudsTemp2 = new int*[Ncpu];
    for (int p = 0; p < Ncpu; p++)
    {
      sommeNoeudsTemp2[p] = new int[nombreFacesTemp[p]];
      facesTemp2[p] = new int*[nombreFacesTemp[p]];
      for (int i = 0; i < nombreFacesTemp[p]; i++)
      {
        facesTemp2[p][i] = new int[4];
      }
    }
    for (int i = nombreElementsFrontiere; i < nombreElementsGlobal; i++)
    {
      if ((i - nombreElementsFrontiere) % frequenceImpression == 0) { cout << "    " << (100 * (i - nombreElementsFrontiere) / (nombreElementsGlobal - nombreElementsFrontiere)) << "% ... " << endl; }
      int numCPU(elementsGlobal[i]->getCPU());
      elementsGlobal[i]->construitFacesSimplifie(iMaxFaces[numCPU], facesTemp2[numCPU], sommeNoeudsTemp2[numCPU]);
    }
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    cout << "    OK en " << t1 << " secondes" << endl;

    //5)Recherche des elements fantomes qui ont une face commune avec un element interne (communicants)
    //-------------------------------------------------------------------------------------------------
    cout << "  5/Recherche elements fantomes ..." << endl;
    tTemp = clock();
    frequenceImpression = max(nombreElementsGlobal / 10, 1);
    int *nombreFacesCommunicantesCPU = new int[Ncpu];
    for (int i = 0; i < Ncpu; i++) { nombreFacesCommunicantesCPU[i] = 0; }
    for (int i = 0; i < nombreElementsGlobal; i++)
    {
      if (i%frequenceImpression == 0) { cout << "    " << (100 * i / nombreElementsGlobal) << "% ... " << endl; }
      if (elementsGlobal[i]->getNombreAutresCPU() != 0)
      {
        vector<int> CPUAEnlever;
        for (int p = 0; p < elementsGlobal[i]->getNombreAutresCPU(); p++)
        {
          int numCPU(elementsGlobal[i]->getAutreCPU(p));
          //verification si l element est communicant via une face
          //******************************************************
          int nombreNoeuds(elementsGlobal[i]->getNombreNoeuds());
          bool *isNoeudInterne = new bool[nombreNoeuds];
          for (int n = 0; n < nombreNoeuds; n++)
          {
            isNoeudInterne[n] = false;
            noeudCourant = elementsGlobal[i]->getNumNoeud(n);
            for (int j = 0; j < nombreNoeudsInterne[numCPU]; j++)
            {
              if (noeudsCPU[numCPU][j] == noeudCourant)
              {
                isNoeudInterne[n] = true;
                break;
              }
            }
          }
          //Determination du nombre de faces communicantes
          //**********************************************
          //int nombreFacesCommunicantes(elementsGlobal[i]->compteFaceCommunicante(facesTemp[numCPU],sommeNoeudsTemp[numCPU]));
          int nombreFacesCommunicantes(elementsGlobal[i]->compteFaceCommunicante(iMaxFaces[numCPU], facesTemp2[numCPU], sommeNoeudsTemp2[numCPU]));
          if (nombreFacesCommunicantes > 0)
          {
            //L'element est communicant, on l'ajoute ainsi que ces noeuds
            elementsCPU[numCPU].push_back(i);
            for (int n = 0; n < nombreNoeuds; n++)
            {
              if (!isNoeudInterne[n]) { noeudsCPU[numCPU].push_back(elementsGlobal[i]->getNumNoeud(n)); }
            }
            nombreFacesCommunicantesCPU[numCPU] += nombreFacesCommunicantes;
          }
          else
          {
            CPUAEnlever.push_back(p);
          }
          delete[] isNoeudInterne;

        } //Fin CPU
        elementsGlobal[i]->enleveCPUAutres(CPUAEnlever);
      }
    }
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    cout << "    OK en " << t1 << " secondes" << endl;

    for (int p = 0; p < Ncpu; p++)
    {
      for (int i = 0; i < nombreFacesTemp[p]; i++) { delete facesTemp2[p][i]; }
      delete sommeNoeudsTemp2[p]; delete facesTemp2[p];
    }
    delete[]facesTemp2; delete[]sommeNoeudsTemp2;
    delete[]nombreFacesTemp; delete[]iMaxFaces;

    //6) Ecriture des fichiers de maillage pour chaque CPU
    cout << "  6/Ecriture des fichiers de maillage pour chacun des " << Ncpu << " CPU ..." << endl;
    tTemp = clock();
    for (int p = 0; p < Ncpu; p++)
    {
      stringstream flux;
      flux << p;
      string fichierMaillageCPU("./libMaillages/" + m_nomMaillage + "_CPU" + flux.str() + ".msh");
      cout << "    ecriture '" << fichierMaillageCPU.c_str() << "' en cours ...' " << endl;
      ofstream fluxFichier;
      fluxFichier.open(fichierMaillageCPU.c_str());

      fluxFichier << "$MeshFormat" << endl;
      fluxFichier << "2.2 0 8" << endl;
      fluxFichier << "$EndMeshFormat" << endl;
      fluxFichier << "$Nodes" << endl;
      //Reordonne les noeuds
      //sort(noeudsCPU[p].begin(), noeudsCPU[p].end());
      fluxFichier << noeudsCPU[p].size() << endl;
      for (unsigned int i = 0; i < noeudsCPU[p].size(); i++)
      {
        fluxFichier << i + 1 << " " << noeudsGlobal[noeudsCPU[p][i]].getX() << " " << noeudsGlobal[noeudsCPU[p][i]].getY() << " " << noeudsGlobal[noeudsCPU[p][i]].getZ() << endl;
      }
      fluxFichier << "$EndNodes" << endl;
      fluxFichier << "$Elements" << endl;
      fluxFichier << elementsCPU[p].size() << endl;
      for (unsigned int i = 0; i < elementsCPU[p].size(); i++)
      {
        int e(elementsCPU[p][i]);
        fluxFichier << i + 1 << " " << elementsGlobal[e]->getTypeGmsh();
        fluxFichier << " " << 2 + 1 + 1 + elementsGlobal[e]->getNombreAutresCPU();
        fluxFichier << " " << elementsGlobal[e]->getAppartenancePhysique();
        fluxFichier << " " << elementsGlobal[e]->getAppartenanceGeometrique();
        fluxFichier << " " << elementsGlobal[e]->getNombreAutresCPU() + 1;
        fluxFichier << " " << elementsGlobal[e]->getCPU() + 1;
        for (int cpuAutre = 0; cpuAutre < elementsGlobal[e]->getNombreAutresCPU(); cpuAutre++)
        {
          fluxFichier << " " << -(elementsGlobal[e]->getAutreCPU(cpuAutre) + 1);
        }
        for (int n = 0; n < elementsGlobal[e]->getNombreNoeuds(); n++)
        {
          //Renumerotation locale
          noeudCourant = elementsGlobal[e]->getNumNoeud(n);
          for (int j = 0; j < static_cast<int>(noeudsCPU[p].size()); j++)
          {
            if (noeudsCPU[p][j] == noeudCourant)
            {
              fluxFichier << " " << j + 1; break;
            }
          }
        }
        fluxFichier << endl;
      }
      fluxFichier << "$EndElements" << endl;
      //Informations additionelles utiles !!
      fluxFichier << "Info non lue par Gmsh : nombre de faces communicante" << endl;
      fluxFichier << nombreFacesCommunicantesCPU[p] << endl;
      fluxFichier << "Info non lue par Gmsh : nombre de noeuds internes (hors fantomes)" << endl;
      fluxFichier << nombreNoeudsInterne[p] << endl;
      fluxFichier.close();
    }
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    cout << "    OK en " << t1 << " secondes" << endl;

    cout << "... FIN PRETRAITEMENT FICHIERS DE MAILLAGE ";
    tempsTotal = clock() - tempsTotal; t1 = static_cast<float>(tempsTotal) / CLOCKS_PER_SEC;
    cout << " Temps de pretraitement total : " << t1 << " secondes" << endl;

    cout << "------------------------------------------------------" << endl;

    //Desallocations
    for (int i = 0; i < nombreElementsGlobal; i++)
    {
      delete elementsGlobal[i];
    }
    delete[] elementsGlobal;
    delete[] noeudsGlobal;
    delete[] nombreNoeudsInterne;
    delete[] nombreFacesCommunicantesCPU;

  }
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************

void MaillageNonStruct::lectureGeometrieGmsh(vector<ElementNS*>** voisinsNoeuds)
{
  try {
    //1) Ouverture du fichier de maillage
    //-------------------------------------
    cout << "------------------------------------------------------" << endl;
    cout << " A) LECTURE DU FICHIER DE MAILLAGE " + m_fichierMaillage + " EN COURS..." << endl;
    m_fichierMaillage = "./libMaillages/" + m_fichierMaillage;
    ifstream fichierMaillage(m_fichierMaillage.c_str(), ios::in);
    if (!fichierMaillage){ throw ErreurECOGEN("fichier maillage absent :" + m_fichierMaillage, __FILE__, __LINE__); }
    string ligneCourante;
    getline(fichierMaillage, ligneCourante);
    getline(fichierMaillage, ligneCourante);
    getline(fichierMaillage, ligneCourante);
    getline(fichierMaillage, ligneCourante);

    //2) Stockage de la grille de point dans tableau m_noeuds
    //-------------------------------------------------------
    cout << "  1/lecture des noeuds du maillage ...";
    fichierMaillage >> m_nombreNoeuds;
    fichierMaillage.ignore(1, '\n');
    m_noeuds = new Coord[m_nombreNoeuds];
    *voisinsNoeuds = new vector<ElementNS*>[m_nombreNoeuds];
    int inutile(0); double x, y, z;
    for (int i = 0; i < m_nombreNoeuds; i++)
    {
      fichierMaillage >> inutile >> x >> y >> z;
      m_noeuds[i].setXYZ(x, y, z);
    }
    fichierMaillage.ignore(1, '\n');
    getline(fichierMaillage, ligneCourante);
    getline(fichierMaillage, ligneCourante);
    cout << "OK" << endl;

    //3) Recuperation des elements 1D/2D/3D dans le tableau m_elements et comptage
    //----------------------------------------------------------------------------
    cout << "  2/lecture des elements 0D/1D/2D/3D ...";
    fichierMaillage >> m_nombreElements;
    fichierMaillage.ignore(1, '\n');
    //Allocation tableau d elements
    m_elements = new ElementNS*[m_nombreElements];
    //Lecture des elements et attributions proprietes geometriques
    m_nombreElements1D = 0, m_nombreElements2D = 0, m_nombreElements3D = 0;
    int posDebutElement = static_cast<int>(fichierMaillage.tellg()); //reperage debut des elements pour retour rapide
    int noeudG;
    for (int i = 0; i < m_nombreElements; i++) {
      this->lectureElementGmsh(m_noeuds, fichierMaillage, &m_elements[i]);
      if (m_elements[i]->getTypeGmsh() == 15) { m_nombreElements0D++; }
      else if (m_elements[i]->getTypeGmsh() == 1) { m_nombreElements1D++; }
      else if (m_elements[i]->getTypeGmsh() <= 3) { m_nombreElements2D++; }
      else if (m_elements[i]->getTypeGmsh() <= 7) { m_nombreElements3D++; }
      else { throw ErreurECOGEN("Type element du .msh non gere dans ECOGEN", __FILE__, __LINE__); }
      //Attribution element i voisin pour les noeuds concernes (Ordre 2 muiltipentes)
      for (int n = 0; n < m_elements[i]->getNombreNoeuds(); n++) {
        noeudG = m_elements[i]->getNumNoeud(n);
        (*voisinsNoeuds)[noeudG].push_back(m_elements[i]);
      }
    }
    m_nombreElementsInternes = m_nombreElements;

    //for (int n = 0; n < m_nombreNoeuds; n++) {
    //  for (int e = 0; e < (*voisinsNoeuds)[n].size(); e++) {
    //    cout << (*voisinsNoeuds)[n][e]->getIndice() << " ";
    //  }
    //  cout << endl;
    //}

    fichierMaillage.close();
    cout << "OK" << endl;
    cout << endl << "  --------------------------" << endl;
    cout << "    INFORMATIONS MAILLAGE :" << endl;
    cout << "  --------------------------" << endl;
    cout << "    nombre de noeuds de maillage : " << m_nombreNoeuds << endl;
    cout << "    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    if (m_nombreElements0D != 0)
    {
      cout << endl << "    nombre elements 0D : " << m_nombreElements0D << endl;
      cout << "    ~~~~~~~~~~~~~~~~~~~~" << endl;
      if (m_nombrePoints != 0) { cout << "      - nombre points : " << m_nombrePoints << endl; }
    }
    if (m_nombreElements1D != 0)
    {
      cout << endl << "    nombre elements 1D : " << m_nombreElements1D << endl;
      cout << "    ~~~~~~~~~~~~~~~~~~~~" << endl;
      if (m_nombreSegments != 0) { cout << "      - nombre segments : " << m_nombreSegments << endl; }
    }
    if (m_nombreElements2D != 0)
    {
      cout << endl << "    nombre elements 2D : " << m_nombreElements2D << endl;
      cout << "    ~~~~~~~~~~~~~~~~~~~~" << endl;
      if (m_nombreTriangles != 0) { cout << "      - nombre triangles   : " << m_nombreTriangles << endl; }
      if (m_nombreQuadrangles != 0) { cout << "      - nombre quadrangles : " << m_nombreQuadrangles << endl; }
    }
    if (m_nombreElements3D != 0)
    {
      cout << endl << "    nombre elements 3D : " << m_nombreElements3D << endl;
      cout << "    ~~~~~~~~~~~~~~~~~~~~" << endl;
      if (m_nombreTetraedres != 0) { cout << "      - nombre tetraedres   : " << m_nombreTetraedres << endl; }
      if (m_nombrePyramides != 0) { cout << "      - nombre pyramides    : " << m_nombrePyramides << endl; }
      if (m_nombreHexaedres != 0) { cout << "      - nombre hexaedres    : " << m_nombreHexaedres << endl; }
    }
    cout << endl << "... FIN LECTURE DU FICHIER DE MAILLAGE " << endl;
    cout << "------------------------------------------------------" << endl;
  }
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************

void MaillageNonStruct::lectureGeometrieGmshParallele()
{
  int nombreNoeudsTotal(0), nombreElementsTotal(0);
  int nombreCPUVoisins(0), CPUCourant(0);

  try{
    //1) Ouverture du fichier de maillage
    //-------------------------------------
    if (rang == 0)
    {
      cout << "------------------------------------------------------" << endl;
      cout << " A) LECTURE DES FICHIERS DE MAILLAGE " + m_nomMaillage + "_CPUX.msh" + " EN COURS..." << endl;
    }
    stringstream flux;
    flux << rang;
    m_fichierMaillage = "./libMaillages/" + m_nomMaillage + "_CPU" + flux.str() + ".msh";
    ifstream fichierMaillage(m_fichierMaillage.c_str(), ios::in);
    if (!fichierMaillage) { throw ErreurECOGEN("fichier maillage absent :" + m_fichierMaillage, __FILE__, __LINE__); }
    string ligneCourante;
    getline(fichierMaillage, ligneCourante);
    getline(fichierMaillage, ligneCourante);
    getline(fichierMaillage, ligneCourante);
    getline(fichierMaillage, ligneCourante);

    //2) Stockage de la grille de point dans tableau m_noeuds
    //-------------------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    if (rang == 0) { cout << "  1/lecture des noeuds du maillage ..."; }
    fichierMaillage >> m_nombreNoeuds;
    fichierMaillage.ignore(1, '\n');
    m_noeuds = new Coord[m_nombreNoeuds];
    int inutile(0); double x, y, z;
    for (int i = 0; i < m_nombreNoeuds; i++)
    {
      fichierMaillage >> inutile >> x >> y >> z;
      m_noeuds[i].setXYZ(x, y, z);
    }
    fichierMaillage.ignore(1, '\n');
    getline(fichierMaillage, ligneCourante);
    getline(fichierMaillage, ligneCourante);
    if (rang == 0) { cout << "OK" << endl; }

    //3) Recuperation des elements 1D/2D/3D dans le tableau m_elements et comptage
    //----------------------------------------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    if (rang == 0) { cout << "  2/lecture des elements interne 1D/2D/3D ..."; }
    fichierMaillage >> m_nombreElements;
    fichierMaillage.ignore(1, '\n');
    //Allocation tableau d elements
    m_elements = new ElementNS*[m_nombreElements];
    //Lecture des elements et attributions proprietes geometriques
    m_nombreElements1D = 0, m_nombreElements2D = 0, m_nombreElements3D = 0;
    int posDebutElement = static_cast<int>(fichierMaillage.tellg()); //reperage debut des elements pour retour rapide
    for (int i = 0; i < m_nombreElements; i++)
    {
      this->lectureElementGmsh(m_noeuds, fichierMaillage, &m_elements[i]);
      if (m_elements[i]->getCPU() == rang)
      {
        if (m_elements[i]->getTypeGmsh() == 1) { m_nombreElements1D++; }
        else if (m_elements[i]->getTypeGmsh() <= 3) { m_nombreElements2D++; }
        else if (m_elements[i]->getTypeGmsh() <= 7) { m_nombreElements3D++; }
        else { throw ErreurECOGEN("Type element du .msh non gere dans ECOGEN", __FILE__, __LINE__); }
      }
      else { m_nombreElementsFantomes++; }
    }
    fichierMaillage.ignore(1, '\n');
    getline(fichierMaillage, ligneCourante);
    //Lecture informations hors Gmsh
    getline(fichierMaillage, ligneCourante);
    fichierMaillage >> m_nombreFacesParallele;
    fichierMaillage.ignore(1, '\n');
    getline(fichierMaillage, ligneCourante);
    fichierMaillage >> m_nombreNoeudsInternes;
    fichierMaillage.close();
    //Calcul du nombre d'elements propres au CPU
    m_nombreElementsInternes = m_nombreElements - m_nombreElementsFantomes;
    if (rang == 0) { cout << "OK" << endl; }
  }
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************

void MaillageNonStruct::lectureElementGmsh(const Coord *TableauNoeuds, ifstream &fichierMaillage, ElementNS **element)
{
  int numeroElement,nombreTags,typeElement,numeroEntitePhysique,numeroEntiteGeometrique;
  bool creeElement(true);
  fichierMaillage >> numeroElement >> typeElement;
  //Lecture des tags
  fichierMaillage >> nombreTags;
  fichierMaillage >> numeroEntitePhysique;    //Appartenance physique
  fichierMaillage >> numeroEntiteGeometrique; //Appartenance geometrique

  //1)Affectation du nombre de point selon element
  //----------------------------------------------
  switch (typeElement)
  {
    case 1: //segment (deux points)
      *element = new ElementSegment;// cout << "segment trouve" << endl;
      m_nombreSegments++;
      break;
    case 2: //triangle (trois points)
      *element = new ElementTriangle;// cout << "triangle trouve" << endl;
      m_nombreTriangles++;
      break;
    case 3: //Quadrangle (quatre points)
      *element = new ElementQuadrangle;// cout << "quadrangle trouve" << endl;
      m_nombreQuadrangles++;
      break;
    case 4: //Tetraedre (quatre points)
      *element = new ElementTetraedre;// cout << "tetraedre trouve" << endl;
      m_nombreTetraedres++;
      break;
    //case 7: //Pyramide quadrangulaire (cinq points) // Ce type d'element semble ne pas fonctionner avec GMSH, les volumes des elements du maillage semblent poser probleme...
    //  *element = new ElementPyramide;
    //  m_nombrePyramides++;
    //  break;
    case 15: //Point (un point)
      *element = new ElementPoint; // cout << "point trouve" << endl;
      m_nombrePoints++;
      break;
    case 5: //Hexaedre (huit points)
      *element = new ElementHexaedre; // cout << "Hexaedre trouve" << endl;
      m_nombreHexaedres++;
      break;
    case 6: //Prisme (six points)
      *element = new ElementPrisme; // cout << "Hexaedre trouve" << endl;
      m_nombreHexaedres++;
      break;
    default:
      Erreurs::messageErreur("Type d element du fichier .msh inconnu de ECOGEN");
      break;
  } //Fin switch typeElement
 
  //2) Specificite maillages paralleles
  //-----------------------------------
  int nombreCPU(0);
  if (nombreTags > 2)
  {
    fichierMaillage >> nombreCPU; //nombre de partition de maillage auquel appartient l element
    int *numCPU = new int[nombreCPU];
    for (int tag = 0; tag < nombreCPU; tag++){ fichierMaillage >> numCPU[tag]; }
    (*element)->setAppartenanceCPU(numCPU, nombreCPU);
    delete[] numCPU;
  }

  //3) Construction de l'element et de ses proprietes
  //-------------------------------------------------
  int noeudCourant;
  int *numNoeud = new int[(*element)->getNombreNoeuds()];
  Coord *noeud = new Coord[(*element)->getNombreNoeuds()];
  for (int i = 0; i < (*element)->getNombreNoeuds(); i++)
  {
    fichierMaillage >> noeudCourant;
    numNoeud[i] = noeudCourant-1;         //decalage car tableau commencant a zero
    noeud[i] = TableauNoeuds[noeudCourant - 1];
  }
  int indiceElement(numeroElement - 1);
  (*element)->construitElement(numNoeud, noeud, numeroEntitePhysique, numeroEntiteGeometrique, indiceElement);

  delete[] noeud;
  delete[] numNoeud;
  
}

//**************************************************************************
//******************************** ECRITURE ********************************
//**************************************************************************

void MaillageNonStruct::ecritEntetePiece(std::ofstream &fluxFichier, std::vector<Cellule *> *cellulesLvl, int lvl) const
{
  fluxFichier << "    <Piece NumberOfPoints=\"" << m_nombreNoeuds << "\" NumberOfCells=\"" << m_nombreCellulesCalcul - m_nombreCellulesFantomes << "\">" << endl;
}

//****************************************************************************

void MaillageNonStruct::recupereNoeuds(std::vector<double> &jeuDonnees, int lvl) const
{
  for (int noeud = 0; noeud < m_nombreNoeuds; noeud++)
  {
    jeuDonnees.push_back(m_noeuds[noeud].getX());
    jeuDonnees.push_back(m_noeuds[noeud].getY());
    jeuDonnees.push_back(m_noeuds[noeud].getZ());
  }
}

//****************************************************************************

void MaillageNonStruct::recupereConnectivite(std::vector<double> &jeuDonnees, int lvl) const
{
  for (int i = m_nombreFacesLimites; i < m_nombreElementsInternes; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      for (int noeud = 0; noeud < m_elements[i]->getNombreNoeuds(); noeud++)
      {
        jeuDonnees.push_back(m_elements[i]->getNumNoeud(noeud));
      }
    }
  }
}

//****************************************************************************

void MaillageNonStruct::recupereOffsets(std::vector<double> &jeuDonnees, int lvl) const
{
  int offset(0);
  for (int i = m_nombreFacesLimites; i < m_nombreElementsInternes; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      offset += m_elements[i]->getNombreNoeuds();
      jeuDonnees.push_back(offset);
    }
  }
}

//****************************************************************************

void MaillageNonStruct::recupereTypeCell(std::vector<double> &jeuDonnees, int lvl) const
{
  for (int i = m_nombreFacesLimites; i < m_nombreElementsInternes; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      jeuDonnees.push_back(m_elements[i]->getTypeVTK());
    }
  }
}

//****************************************************************************
/*!
* \brief     Recuperation de donnees pour ecriture
* \details   Fonctions permettant de recuperer un esemble de donnees de phase ou de melange,
*                  scalaire ou vectorielle
* \param    cellules         cellules contenant les donnees
* \param    var              numero de la variable demandee (>0 si scalaire, <0 si vectorielle)
* \param    phase            numero de la phase demandee (-1 si melange demande, -2 pour transports, -3 pour xi, -4 pour gradient densite melange)
* \return   jeuDonnees       vecteur de double contenant toutes les donnees de la variables var de la phase phase a la suite
*/
void MaillageNonStruct::recupereDonnees(vector<Cellule *> *cellulesLvl, std::vector<double> &jeuDonnees, const int var, int phase, int lvl) const
{
  int numCell;
  for (int i = m_nombreFacesLimites; i < m_nombreElementsInternes; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      numCell = m_elements[i]->getNumCelluleAssociee();
      if (var > 0) { //On veut recuperer les donnees scalaires
        if (phase >= 0) { jeuDonnees.push_back(cellulesLvl[0][numCell]->getPhase(phase)->renvoieScalaire(var)); } //Donnees de phases
        else if (phase == -1) { jeuDonnees.push_back(cellulesLvl[0][numCell]->getMelange()->renvoieScalaire(var)); }               //Donnees de melange
        else if (phase == -2) { jeuDonnees.push_back(cellulesLvl[0][numCell]->getTransport(var-1).getValeur()); }
        else if (phase == -3) { jeuDonnees.push_back(cellulesLvl[0][numCell]->getXi()); }
        else if (phase == -4) { jeuDonnees.push_back(cellulesLvl[0][numCell]->getGradient()); }
        else { Erreurs::messageErreur("MaillageNonStruct::recupereDonnees : numero de phase inconnu : ", phase); }
      }
      else { //On veut recuperer les donnees vectorielles
        if (phase >= 0) { //Donnees de phases
          jeuDonnees.push_back(cellulesLvl[0][numCell]->getPhase(phase)->renvoieVecteur(-var)->getX());
          jeuDonnees.push_back(cellulesLvl[0][numCell]->getPhase(phase)->renvoieVecteur(-var)->getY());
          jeuDonnees.push_back(cellulesLvl[0][numCell]->getPhase(phase)->renvoieVecteur(-var)->getZ());
        }
        else if(phase == -1){  //Donnees de melange
          jeuDonnees.push_back(cellulesLvl[0][numCell]->getMelange()->renvoieVecteur(-var)->getX());
          jeuDonnees.push_back(cellulesLvl[0][numCell]->getMelange()->renvoieVecteur(-var)->getY());
          jeuDonnees.push_back(cellulesLvl[0][numCell]->getMelange()->renvoieVecteur(-var)->getZ());
        }
        else { Erreurs::messageErreur("MaillageNonStruct::recupereDonnees : numero de phase inconnu : ", phase); }
      } //Fin vecteur
    }
  }
}

//***********************************************************************

void MaillageNonStruct::rechercheElementsArrieres(ElementNS *element, FaceNS *face, BordDeMaille *bord, vector<ElementNS *> voisins, Cellule **cellules) const
{
  //int indice(element->getIndice() - m_nombreFacesLimites);
  //Coord vG = element->vecteur(face);
  //ElementNS *eT; //Element a tester
  ////1) Recherche du premier element arriere gauche BG1M
  ////---------------------------------------------------
  //double cos(0.);
  //for (unsigned int v = 0; v < voisins.size(); v++) {
  //  eT = voisins[v];
  //  Coord vT = eT->vecteur(element);
  //  double cosTemp(Coord::cos(vT, vG));
  //  if (cosTemp >= cos) {
  //    cos = cosTemp;
  //    if (eT->getIndice() < m_nombreFacesLimites) { //Si c est une limite on remet BG1M a NULL
  //      bord->setB(BG1M, 0);
  //    }
  //    else { //Sinon on met a jour avec la maille la plus loin
  //      Cellule *c(cellules[eT->getNumCelluleAssociee()]);
  //      bord->setB(BG1M, c);
  //    }
  //  }
  //}
  ////2) Recherche du second element arriere gauche BG2M
  ////---------------------------------------------------
  //cos = 0.;
  //if (bord->getB(BG1M) != 0) {  //Cas non bord
  //  Element *e1(bord->getB(BG1M)->getElement());
  //  Coord v1 = e1->vecteur(element);
  //  Coord sin1(Coord::sin(v1, vG));
  //  if (fabs(sin1.norme()) <= 1e-8) {
  //    bord->setB(BG2M, bord->getB(BG1M)); // Si le cos est nul, meme cellule pour B1 et B2
  //  }
  //  else {
  //    for (unsigned int v = 0; v < voisins.size(); v++) {
  //      eT = voisins[v];
  //      if (eT == e1) continue; // voisin suivant si voisin deja utilise 
  //      Coord vT = eT->vecteur(element);
  //      Coord sinTemp(Coord::sin(vT, vG));

  //      if (sinTemp.scalaire(sin1) <= 0.) {
  //        double cosTemp(Coord::cos(vT, vG));
  //        if (cosTemp >= cos) {
  //          cos = cosTemp;
  //          if (eT->getIndice() < m_nombreFacesLimites) { //Si c est une limite on remet BG1M a NULL
  //            bord->setB(BG2M, 0);
  //          }
  //          else {  //Sinon on met a jour avec la 2 eme maille la plus loin
  //            Cellule *c(cellules[eT->getNumCelluleAssociee()]);
  //            bord->setB(BG2M, c);
  //          }
  //        }
  //      } //fin sin*sin <0
  //    }  //fin boucle voisins
  //  } //fin if 
  //} //Fin if bord

  //  //Determination des ponderations arrieres
  //Coord MB1; if (bord->getB(BG1M) != 0) { MB1.creeVecteur(face->getPos(), bord->getB(BG1M)->getPosition()); }
  //Coord MB2; if (bord->getB(BG2M) != 0) { MB2.creeVecteur(face->getPos(), bord->getB(BG2M)->getPosition()); }
  //Coord MB; MB.creeVecteur(face->getPos(), element->getPosition());
  //double a, b, beta1(1.), beta2(0.);
  //a = (MB1.vectoriel(MB)).norme() / MB.norme();
  //b = (MB2.vectoriel(MB)).norme() / MB.norme();
  //if (fabs(a + b) > 1e-8) { beta1 = b / (a + b); beta2 = 1. - beta1; }
  //bord->setBeta(betaG1M, beta1);
  //bord->setBeta(betaG2M, beta2);
  ////Calcul de la distance du point M (bord de maille) au point H (barycentre)
  //if (bord->getB(BG1M) != 0) {
  //  cos = Coord::cos(MB, MB1);
  //  double d, e, c;
  //  d = cos*MB1.norme();
  //  Coord B1B2; B1B2 = MB1 - MB2;
  //  e = B1B2.norme()*beta1;
  //  c = sqrt(e*e - a*a);
  //  d += c;
  //  bord->setDistanceH(distanceHGM, d);
  //}
  //cout << eG->getIndice() << " : " << endl;
  //if ((*bord)[i]->getB(BG1M) != 0) cout << "BG1M = " << (*bord)[i]->getB(BG1M)->getElement()->getIndice() << " " << (*bord)[i]->getDistanceH(distanceHGM) << endl;
  //if ((*bord)[i]->getB(BG2M) != 0) cout << "BG2M = " << (*bord)[i]->getB(BG2M)->getElement()->getIndice() << " " << (*bord)[i]->getDistanceH(distanceHGM) << endl;
}

//***********************************************************************

void MaillageNonStruct::rechercheElementsAvants(ElementNS *element, FaceNS *face, BordDeMaille *bord, vector<ElementNS *> voisins, Cellule **cellules) const
{
  int indice(element->getIndice() - m_nombreFacesLimites);
  Coord vG = element->vecteur(face);
  //ElementNS *eT; //Element a tester
  ////1) Recherche du premier element arriere gauche BG1P
  ////---------------------------------------------------
  //double cos(0.);
  //for (unsigned int v = 0; v < voisins.size(); v++) {
  //  eT = voisins[v];
  //  Coord vT = eT->vecteur(element);
  //  double cosTemp(Coord::cos(vT, vG));
  //  if (cosTemp <= cos) {
  //    cos = cosTemp;
  //    if (eT->getIndice() < m_nombreFacesLimites) { //Si c est une limite on remet BG1P a NULL
  //      bord->setB(BG1P, 0);
  //    }
  //    else { //Sinon on met a jour avec la maille la plus loin
  //      Cellule *c(cellules[eT->getNumCelluleAssociee()]);
  //      bord->setB(BG1P, c);
  //    }
  //  }
  //}
  ////2) Recherche du second element arriere gauche BG2P
  ////---------------------------------------------------
  //cos = 0.;
  //if (bord->getB(BG1P) != 0) {  //Cas non bord
  //  Element *e1(bord->getB(BG1P)->getElement());
  //  Coord v1 = e1->vecteur(element);
  //  Coord sin1(Coord::sin(v1, vG));
  //  if (fabs(sin1.norme()) <= 1e-8) {
  //    bord->setB(BG2P, bord->getB(BG1P)); // Si le cos est nul, meme cellule pour B1 et B2
  //  }
  //  else {
  //    for (unsigned int v = 0; v < voisins.size(); v++) {
  //      eT = voisins[v];
  //      if (eT == e1) continue; // voisin suivant si voisin deja utilise 
  //      Coord vT = eT->vecteur(element);
  //      Coord sinTemp(Coord::sin(vT, vG));

  //      if (sinTemp.scalaire(sin1) <= 0.) {
  //        double cosTemp(Coord::cos(vT, vG));
  //        if (cosTemp <= cos) {
  //          cos = cosTemp;
  //          if (eT->getIndice() < m_nombreFacesLimites) { //Si c est une limite on remet BG2P a NULL
  //            bord->setB(BG2P, 0);
  //          }
  //          else {  //Sinon on met a jour avec la 2 eme maille la plus loin
  //            Cellule *c(cellules[eT->getNumCelluleAssociee()]);
  //            bord->setB(BG2P, c);
  //          }
  //        }
  //      } //fin sin*sin <0
  //    }  //fin boucle voisins
  //  } //fin if 
  //} //Fin if bord

  //  //Determination des ponderations arrieres
  //Coord MB1; if (bord->getB(BG1P) != 0) { MB1.creeVecteur(face->getPos(), bord->getB(BG1P)->getPosition()); }
  //Coord MB2; if (bord->getB(BG2P) != 0) { MB2.creeVecteur(face->getPos(), bord->getB(BG1P)->getPosition()); }
  //Coord MB; MB.creeVecteur(face->getPos(), element->getPosition());
  //double a, b, beta1(1.), beta2(0.);
  //a = (MB1.vectoriel(MB)).norme() / MB.norme();
  //b = (MB2.vectoriel(MB)).norme() / MB.norme();
  //if (fabs(a + b) > 1e-8) { beta1 = b / (a + b); beta2 = 1. - beta1; }
  //bord->setBeta(betaG1P, beta1);
  //bord->setBeta(betaG2P, beta2);
  ////Calcul de la distance du point M (bord de maille) au point H (barycentre)
  //if (bord->getB(BG1P) != 0) {
  //  cos = Coord::cos(MB, -1.*MB1);
  //  double d, e, c;
  //  d = cos*MB1.norme();
  //  Coord B1B2; B1B2 = MB1 - MB2;
  //  e = B1B2.norme()*beta1;
  //  c = sqrt(e*e - a*a);
  //  d += c;
  //  bord->setDistanceH(distanceHGP, d);
  //}
  //cout << eG->getIndice() << " : " << endl;
  //if ((*bord)[i]->getB(BG1P) != 0) cout << "BG1P = " << (*bord)[i]->getB(BG1P)->getElement()->getIndice() << " " << (*bord)[i]->getDistanceH(distanceHGP) << endl;
  //if ((*bord)[i]->getB(BG2P) != 0) cout << "BG2P = " << (*bord)[i]->getB(BG2P)->getElement()->getIndice() << " " << (*bord)[i]->getDistanceH(distanceHGP) << endl;
}

//************************REPRISE DEPUIS FICHIERS RESULTATS***************************************

void MaillageNonStruct::repriseCalcul(Cellule **cellules, const int &nombrePhases, const int &nombreTransports, const int &repriseFichier, bool ecritVTK, bool ecritXML, bool ecritBinaire, bool cree)
{
  m_numFichier = repriseFichier; //On attribut le numero de fichier

  if (ecritVTK&&ecritXML&&ecritBinaire) {
    /*if (Ncpu == 1) { this->ecritSolutionXMLBinaire(cellules, nombrePhases); }
    else { this->ecritSolutionXMLParallelBinaire(cellules, nombrePhases); }*/
    Erreurs::messageErreur("Format inconnu sur maillage non structure pour la reprise depuis fichier");
  }
  else if (ecritVTK&&ecritXML&&!ecritBinaire) {
    //if (Ncpu == 1) { this->ecritSolutionXML(cellules, nombrePhases); }
    //else { this->ecritSolutionXMLParallel(cellules, nombrePhases); }
    this->repriseXML(cellules, nombrePhases);
  }
  else if (ecritVTK&&!ecritXML&&ecritBinaire) {
    //this->ecritSolutionVTKBinaire(cellules, nombrePhases);
    Erreurs::messageErreur("Format inconnu sur maillage non structure pour la reprise depuis fichier");
  }
  else if (ecritVTK&&!ecritXML&&!ecritBinaire) {
    //this->ecritSolutionVTK(cellules, nombrePhases);
    Erreurs::messageErreur("Format inconnu sur maillage non structure pour la reprise depuis fichier");
  }
  else { Erreurs::messageErreur("Format inconnu sur maillage non structure pour la reprise depuis fichier"); }

  //On complete les vecteurs ?

  MPI_Barrier(MPI_COMM_WORLD);
  if (rang == 0) cout << "Initialisation OK depuis fichier : " << m_numFichier << endl;
  m_numFichier++;
}

//***********************************************************************

void MaillageNonStruct::repriseXML(Cellule **cellules, const int &nombrePhases) const
{
  int nombreCelluleInterne(m_nombreCellulesCalcul - m_nombreCellulesFantomes);

  stringstream nomFichier;
  nomFichier << "./VTK/result_" << "_CPU_" << rang << "_" << m_numFichier << ".vtu";
  XMLDocument xmlResult;
  XMLError erreur(xmlResult.LoadFile(nomFichier.str().c_str())); //Le fichier est parse ici
  if (erreur != XML_SUCCESS) throw ErreurXML(nomFichier.str(), __FILE__, __LINE__);
  
  //Entree dans l'arborescence
  XMLNode *vtkFILE = xmlResult.FirstChildElement("VTKFile");
  if (vtkFILE == NULL) throw ErreurXMLRacine("VTKFile", nomFichier.str(), __FILE__, __LINE__);
  XMLElement *element0, *element1;
  element0 = vtkFILE->FirstChildElement("UnstructuredGrid");
  if (element0 == NULL) throw ErreurXMLElement("UnstructuredGrid", nomFichier.str(), __FILE__, __LINE__);
  element1 = element0->FirstChildElement("Piece");
  if (element1 == NULL) throw ErreurXMLElement("Piece", nomFichier.str(), __FILE__, __LINE__);
  XMLElement *donnees;
  donnees = element1->FirstChildElement("CellData");
  if (donnees == NULL) throw ErreurXMLElement("CellData", nomFichier.str(), __FILE__, __LINE__);

  int numCell;
  XMLElement *dataArray;
  dataArray = donnees->FirstChildElement("DataArray");
  stringstream fluxValeurs;
  double valeur;
  for (int phase = 0; phase < nombrePhases; phase++)
  {
    //Lecture des variables scalaires
    for (int var = 1; var <= cellules[0]->getPhase(phase)->getNombreScalaires(); var++)
    {
      fluxValeurs.flush();
      fluxValeurs << dataArray->GetText(); //On mets les valeurs dans un flux
      for (int i = m_nombreFacesLimites; i < m_nombreElementsInternes; i++)
      {
        if (!m_elements[i]->isFantome())
        {
          numCell = m_elements[i]->getNumCelluleAssociee();
          fluxValeurs >> valeur;
          cellules[numCell]->getPhase(phase)->setScalaire(var,valeur);
        }
      }
      dataArray = dataArray->NextSiblingElement("DataArray");
    }

    //La temperature ne necessite pas d etre recuperee
    dataArray = dataArray->NextSiblingElement("DataArray");
    
    //Ecriture des variables vectorielles
    for (int var = 1; var <= cellules[0]->getPhase(phase)->getNombreVecteurs(); var++)
    {
      fluxValeurs << dataArray->GetText(); //On mets les valeurs dans un flux
      for (int i = m_nombreFacesLimites; i < m_nombreElementsInternes; i++)
      {
        if (!m_elements[i]->isFantome())
        {
          numCell = m_elements[i]->getNumCelluleAssociee();
          fluxValeurs >> valeur;
          cellules[numCell]->getPhase(phase)->renvoieVecteur(var)->setX(valeur); //FP//ERR// pas bon ici
          fluxValeurs >> valeur;
          cellules[numCell]->getPhase(phase)->renvoieVecteur(var)->setY(valeur);
          fluxValeurs >> valeur;
          cellules[numCell]->getPhase(phase)->renvoieVecteur(var)->setZ(valeur);
        }
      }
      dataArray = dataArray->NextSiblingElement("DataArray");
    }
  } //Fin phases
}

//***********************************************************************