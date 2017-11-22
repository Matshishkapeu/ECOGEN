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

#include "ElementTetraedre.h"

using namespace std;

const int ElementTetraedre::TYPEGMSH = 4;
const int ElementTetraedre::NOMBRENOEUDS = 4;
const int ElementTetraedre::NOMBREFACES = 4; /* ici il s'agit de triangles*/
const int ElementTetraedre::TYPEVTK = 10;

//***********************************************************************

ElementTetraedre::ElementTetraedre() :
ElementNS(TYPEGMSH, NOMBRENOEUDS, NOMBREFACES, TYPEVTK)
{}

//***********************************************************************

ElementTetraedre::~ElementTetraedre(){}

//***********************************************************************

void ElementTetraedre::calculVolume(const Coord *noeuds)
{
  Coord v1(noeuds[1] - noeuds[0]), v2(noeuds[2] - noeuds[0]), v3(noeuds[3] - noeuds[0]);
  m_volume = fabs(Coord::determinant(v1, v2, v3)) / 6.; //volume du tetradre
}

//***********************************************************************

void ElementTetraedre::calculLCFL(const Coord *noeuds)
{
  Coord vec; m_lCFL = 1e10;
  vec = ((noeuds[0] + noeuds[1] + noeuds[2]) / 3.) - m_position;
  m_lCFL = min(m_lCFL, vec.norme());
  vec = ((noeuds[1] + noeuds[2] + noeuds[3]) / 3.) - m_position;
  m_lCFL = min(m_lCFL, vec.norme());
  vec = ((noeuds[2] + noeuds[3] + noeuds[0]) / 3.) - m_position;
  m_lCFL = min(m_lCFL, vec.norme());
  vec = ((noeuds[3] + noeuds[0] + noeuds[1]) / 3.) - m_position;
  m_lCFL = min(m_lCFL, vec.norme());
}

//***********************************************************************
// Nouvelle version beaucoup plus efficace avec recherche dans tableau temporaire
void ElementTetraedre::construitFaces(const Coord *noeuds, FaceNS **faces, int &iMax, int** facesTemp, int* sommeNoeudsTemp)
{
  //4 faces a traiter de type triangle
  int indiceFaceExiste(-1);
  int noeudAutre;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; facesTemp[iMax][2] = m_numNoeuds[2]; noeudAutre = 3; break;
      case 1: facesTemp[iMax][0] = m_numNoeuds[1]; facesTemp[iMax][1] = m_numNoeuds[2]; facesTemp[iMax][2] = m_numNoeuds[3]; noeudAutre = 0; break;      
      case 2: facesTemp[iMax][0] = m_numNoeuds[2]; facesTemp[iMax][1] = m_numNoeuds[3]; facesTemp[iMax][2] = m_numNoeuds[0]; noeudAutre = 1; break;      
      case 3: facesTemp[iMax][0] = m_numNoeuds[3]; facesTemp[iMax][1] = m_numNoeuds[0]; facesTemp[iMax][2] = m_numNoeuds[1]; noeudAutre = 2; break;      
    }
    sommeNoeudsTemp[iMax] = facesTemp[iMax][0] + facesTemp[iMax][1] + facesTemp[iMax][2];
    std::sort(facesTemp[iMax],facesTemp[iMax]+3);  //Tri des noeuds
    //Existance face ?
    indiceFaceExiste = FaceNS::rechercheFace(facesTemp[iMax],sommeNoeudsTemp[iMax],facesTemp,sommeNoeudsTemp,3,iMax);
    //Creation face ou rattachement
    if (indiceFaceExiste==-1)
    {
      faces[iMax] = new FaceTriangle(facesTemp[iMax][0], facesTemp[iMax][1], facesTemp[iMax][2], 0); //pas besoin du tri ici
      faces[iMax]->construitFace(noeuds, m_numNoeuds[noeudAutre], this);
      iMax++;
    }
    else
    {
      faces[indiceFaceExiste]->ajouteElementVoisin(this);
    }
  }
}

//***********************************************************************
// Nouvelle version beaucoup plus efficace avec recherche dans tableau temporaire
void ElementTetraedre::construitFacesSimplifie(int &iMax, int** facesTemp, int* sommeNoeudsTemp)
{
  //4 faces a traiter de type triangle
  int indiceFaceExiste(-1);
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; facesTemp[iMax][2] = m_numNoeuds[2]; break;
      case 1: facesTemp[iMax][0] = m_numNoeuds[1]; facesTemp[iMax][1] = m_numNoeuds[2]; facesTemp[iMax][2] = m_numNoeuds[3]; break;      
      case 2: facesTemp[iMax][0] = m_numNoeuds[2]; facesTemp[iMax][1] = m_numNoeuds[3]; facesTemp[iMax][2] = m_numNoeuds[0]; break;      
      case 3: facesTemp[iMax][0] = m_numNoeuds[3]; facesTemp[iMax][1] = m_numNoeuds[0]; facesTemp[iMax][2] = m_numNoeuds[1]; break;      
    }
    sommeNoeudsTemp[iMax] = facesTemp[iMax][0] + facesTemp[iMax][1] + facesTemp[iMax][2];
    std::sort(facesTemp[iMax],facesTemp[iMax]+3);  //Tri des noeuds
    //Existance face ?
    indiceFaceExiste = FaceNS::rechercheFace(facesTemp[iMax],sommeNoeudsTemp[iMax],facesTemp,sommeNoeudsTemp,3,iMax);
    //Creation face ou rattachement
    if (indiceFaceExiste==-1)
    {
      iMax++;
    }
  }
}

//***********************************************************************

void ElementTetraedre::attributFaceCommunicante(const Coord *noeuds, FaceNS **faces, const int &indiceMaxFaces, const int &nombreNoeudsInternes)
{
  int indiceFaceExiste(0);
  //Verification face 1 :
  if (m_numNoeuds[0] < nombreNoeudsInternes && m_numNoeuds[1] < nombreNoeudsInternes && m_numNoeuds[2] < nombreNoeudsInternes)
  {
    FaceTriangle face(m_numNoeuds[0], m_numNoeuds[1], m_numNoeuds[2]);
    if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
    {
      faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indiceFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 2 :
  if (m_numNoeuds[1] < nombreNoeudsInternes && m_numNoeuds[2] < nombreNoeudsInternes && m_numNoeuds[3] < nombreNoeudsInternes)
  {
    FaceTriangle face(m_numNoeuds[1], m_numNoeuds[2], m_numNoeuds[3]);
    if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
    {
      faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indiceFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 3 :
  if (m_numNoeuds[2] < nombreNoeudsInternes && m_numNoeuds[3] < nombreNoeudsInternes && m_numNoeuds[0] < nombreNoeudsInternes)
  {
    FaceTriangle face(m_numNoeuds[2], m_numNoeuds[3], m_numNoeuds[0]);
    if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
    {
      faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indiceFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 4 :
  if (m_numNoeuds[3] < nombreNoeudsInternes && m_numNoeuds[0] < nombreNoeudsInternes && m_numNoeuds[1] < nombreNoeudsInternes)
  {
    FaceTriangle face(m_numNoeuds[3], m_numNoeuds[0], m_numNoeuds[1]);
    if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
    {
      faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indiceFaceExiste]->setEstComm(true);
    }
  }
}

//***********************************************************************

int ElementTetraedre::compteFaceCommunicante(vector<int*> &facesTemp, vector<int> &sommeNoeudsTemp)
{
  //4 faces a traiter de type triangle
  int indiceFaceExiste(-1), nombreFacesCommunicante(0);
  int face[3], sommeNoeuds;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; face[2] = m_numNoeuds[2]; break;
      case 1: face[0] = m_numNoeuds[1]; face[1] = m_numNoeuds[2]; face[2] = m_numNoeuds[3]; break;
      case 2: face[0] = m_numNoeuds[2]; face[1] = m_numNoeuds[3]; face[2] = m_numNoeuds[0]; break;     
      case 3: face[0] = m_numNoeuds[3]; face[1] = m_numNoeuds[0]; face[2] = m_numNoeuds[1]; break;     
    }
    int iMax = sommeNoeudsTemp.size();
    sommeNoeuds = face[0]+face[1]+face[2];
    std::sort(face, face+3);
    //Recherche existance faces
    indiceFaceExiste = FaceNS::rechercheFace(face,sommeNoeuds,facesTemp,sommeNoeudsTemp,3,iMax);
    if (indiceFaceExiste!=-1)
    {
      nombreFacesCommunicante++;
    }
  }
  return nombreFacesCommunicante;
}

//***********************************************************************
//Nouvelle version plus efficace
int ElementTetraedre::compteFaceCommunicante(int &iMax, int **facesTemp, int *sommeNoeudsTemp)
{
  //4 faces a traiter de type triangle
  int indiceFaceExiste(-1), nombreFacesCommunicante(0);
  int face[3], sommeNoeuds;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; face[2] = m_numNoeuds[2]; break;
      case 1: face[0] = m_numNoeuds[1]; face[1] = m_numNoeuds[2]; face[2] = m_numNoeuds[3]; break;
      case 2: face[0] = m_numNoeuds[2]; face[1] = m_numNoeuds[3]; face[2] = m_numNoeuds[0]; break;     
      case 3: face[0] = m_numNoeuds[3]; face[1] = m_numNoeuds[0]; face[2] = m_numNoeuds[1]; break;     
    }
    sommeNoeuds = face[0]+face[1]+face[2];
    std::sort(face, face+3);
    //Recherche existance faces
    indiceFaceExiste = FaceNS::rechercheFace(face,sommeNoeuds,facesTemp,sommeNoeudsTemp,3,iMax);
    if (indiceFaceExiste!=-1)
    {
      nombreFacesCommunicante++;
    }
  }
  return nombreFacesCommunicante;
}

//***********************************************************************