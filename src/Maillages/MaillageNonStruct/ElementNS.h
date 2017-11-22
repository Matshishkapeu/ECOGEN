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

#ifndef ELEMENTNS_H
#define ELEMENTNS_H

#include "../Element.h"

class FaceNS; //Predeclaration de la classe Face pour pouvoir inclure Face.h
#include "FaceNS.h"

class ElementNS : public Element
{
public:
  ElementNS();
  ElementNS(const int &typeGmsh, const int &nombreNoeuds, const int &nombreFaces, const int &typeVTK);
  virtual ~ElementNS();

  void construitElement(const int *numNoeuds, const Coord *noeuds, const int numeroEntitePhysique, const int numeroEntiteGeometrique, int &indiceElement); /*Calcul des proprietes de l element*/
  void construitElementParallele(const Coord *noeuds); /*Calcul des proprietes de l element*/
  virtual void attributFaceLimite(const Coord *noeuds, FaceNS **faces, const int &indiceMaxFaces){ Erreurs::messageErreur("attributFaceLimite non prevu pour le type d element demande"); };
  virtual void attributFaceCommunicante(const Coord *noeuds, FaceNS **faces, const int &indiceMaxFaces, const int &nombreNoeudsInternes){ Erreurs::messageErreur("attributFaceCommunicante non prevu pour le type d element demande"); };
  virtual void construitFaces(const Coord *noeuds, FaceNS **faces, int &indiceMaxFaces, int** facesTemp, int* sommeNoeudsTemp) { Erreurs::messageErreur("construitFaces non prevu pour le type d element demande"); }; //Pour tests
  virtual void construitFacesSimplifie(int &iMax, int** facesTemp, int* sommeNoeudsTemp){ Erreurs::messageErreur("construitFacesSimplifie non prevu pour le type d element demande"); };
  virtual int compteFaceCommunicante(std::vector<int*> &faces, std::vector<int> &sommeNoeudsTemp){ Erreurs::messageErreur("compteFaceCommunicante non prevu pour le type d element demande"); return 0; };
  virtual int compteFaceCommunicante(int &iMax, int **faces, int *sommeNoeudsTemp){ Erreurs::messageErreur("compteFaceCommunicante non prevu pour le type d element demande"); return 0; };
  /*  void enleveCPUAutres(const int &numCPU);*/
  void enleveCPUAutres(std::vector<int> &numCPU);

  //Accesseurs
  void setIndice(int &indice);
  void setAppartenancePhysique(int &appartenancePhysique);
  void setNumNoeud(int *numNoeuds);
  void setNumNoeud(int &noeud, int &numNoeud);
  void setIsFantome(bool isFantome);
  void setIsCommunicant(bool isCommunicant);
  void setAppartenanceCPU(const int *numCPU, const int &nombreCPU);

  virtual int getIndice() const;
  int getNombreNoeuds() const;
  int getNombreFaces() const;
  int getTypeGmsh() const;
  int getTypeVTK() const;
  int getNumNoeud(int &noeud) const;
  int getAppartenancePhysique() const;
  int getAppartenanceGeometrique() const;
  int getCPU() const;
  int getNombreAutresCPU() const;
  int getAutreCPU(const int &autreCPU) const;
  void afficheInfos() const;

  bool isFantome() const;
  bool isCommunicant() const;

protected:
  virtual void calculVolume(const Coord *noeuds){};
  virtual void calculLCFL(const Coord *noeuds){};

  int m_indice;

  int m_typeGmsh;
  int m_typeVTK;
  int m_nombreNoeuds;
  int m_nombreFaces;
  int m_appartenancePhysique;
  int m_appartenanceGeometrique;
  bool m_isFantome;
  bool m_isCommunicant;
  int m_CPU;              /*Numero du CPU sur lequel l'element est present physiquement*/
  int m_nombreautresCPU;
  int *m_autresCPU;       /*Numero des CPU sur lesquels l'element est present en tant que fantome*/
  int *m_numNoeuds;       /*Correspondance avec le tableau de noeuds du maillage*/

};

#endif // ELEMENTNS_H