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

#include "ElementPrisme.h"

using namespace std;

const int ElementPrisme::TYPEGMSH = 6;
const int ElementPrisme::NOMBRENOEUDS = 6;
const int ElementPrisme::NOMBREFACES = 5; /* ici il s'agit de 3 quadrangles et 2 triangles*/
const int ElementPrisme::TYPEVTK = 13;

//***********************************************************************

ElementPrisme::ElementPrisme() :
ElementNS(TYPEGMSH, NOMBRENOEUDS, NOMBREFACES, TYPEVTK)
{}

//***********************************************************************

ElementPrisme::~ElementPrisme(){}

//***********************************************************************

void ElementPrisme::calculVolume(const Coord *noeuds)
{
  //On va calculer le volume des 3 tetraedres inclus dans le prisme
  Coord v1, v2, v3;
  v1.creeVecteur(noeuds[0], noeuds[1]); v2.creeVecteur(noeuds[0], noeuds[2]); v3.creeVecteur(noeuds[0], noeuds[3]);
  double volumeT1 = fabs(Coord::determinant(v1, v2, v3)) / 6.; //volume du tetradre
  v1.creeVecteur(noeuds[3], noeuds[4]); v2.creeVecteur(noeuds[3], noeuds[5]); v3.creeVecteur(noeuds[3], noeuds[2]);
  double volumeT2 = fabs(Coord::determinant(v1, v2, v3)) / 6.; //volume du tetradre
  v1.creeVecteur(noeuds[3], noeuds[4]); v2.creeVecteur(noeuds[3], noeuds[2]); v3.creeVecteur(noeuds[3], noeuds[1]);
  double volumeT3 = fabs(Coord::determinant(v1, v2, v3)) / 6.; //volume du tetradre
  m_volume = volumeT1 + volumeT2 + volumeT3; //volume du prisme
}

//***********************************************************************

void ElementPrisme::calculLCFL(const Coord *noeuds)
{
  Coord vec; m_lCFL = 1e10;
  vec = ((noeuds[0] + noeuds[1] + noeuds[4] + noeuds[3]) / 4.) - m_position;
  m_lCFL = min(m_lCFL, vec.norme());
  vec = ((noeuds[0] + noeuds[2] + noeuds[5] + noeuds[3]) / 4.) - m_position;
  m_lCFL = min(m_lCFL, vec.norme());
  vec = ((noeuds[1] + noeuds[2] + noeuds[5] + noeuds[4]) / 4.) - m_position;
  m_lCFL = min(m_lCFL, vec.norme());
  vec = ((noeuds[0] + noeuds[1] + noeuds[2]) / 3.) - m_position;
  m_lCFL = min(m_lCFL, vec.norme());
  vec = ((noeuds[3] + noeuds[4] + noeuds[5]) / 3.) - m_position;
  m_lCFL = min(m_lCFL, vec.norme());
}

//***********************************************************************
// Nouvelle version beaucoup plus efficace avec recherche dans tableau temporaire
void ElementPrisme::construitFaces(const Coord *noeuds, FaceNS **faces, int &iMax, int** facesTemp, int* sommeNoeudsTemp)
{
  //3 faces a traiter de type quadrangle et 2 faces de type triangle
  int indiceFaceExiste(-1);
  int noeudAutre;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; facesTemp[iMax][2] = m_numNoeuds[4]; facesTemp[iMax][3] = m_numNoeuds[3]; noeudAutre = 2; break;
      case 1: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[2]; facesTemp[iMax][2] = m_numNoeuds[5]; facesTemp[iMax][3] = m_numNoeuds[3]; noeudAutre = 1; break;      
      case 2: facesTemp[iMax][0] = m_numNoeuds[1]; facesTemp[iMax][1] = m_numNoeuds[2]; facesTemp[iMax][2] = m_numNoeuds[5]; facesTemp[iMax][3] = m_numNoeuds[4]; noeudAutre = 0; break;      
      case 3: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; facesTemp[iMax][2] = m_numNoeuds[2]; noeudAutre = 3; break;      
      case 4: facesTemp[iMax][0] = m_numNoeuds[3]; facesTemp[iMax][1] = m_numNoeuds[4]; facesTemp[iMax][2] = m_numNoeuds[5]; noeudAutre = 0; break;      
    }
    if (i < 3) //Faces Quadrangles
    {
      sommeNoeudsTemp[iMax] = facesTemp[iMax][0] + facesTemp[iMax][1] + facesTemp[iMax][2] + facesTemp[iMax][3];
      std::sort(facesTemp[iMax],facesTemp[iMax]+4);  //Tri des noeuds
      //Existance face ?
      indiceFaceExiste = FaceNS::rechercheFace(facesTemp[iMax],sommeNoeudsTemp[iMax],facesTemp,sommeNoeudsTemp,4,iMax);
      //Creation face ou rattachement
      if (indiceFaceExiste==-1) // on rempli simultanement le tableau faces et le tableau facesTemp
      {
        faces[iMax] = new FaceQuadrangle(facesTemp[iMax][0], facesTemp[iMax][1], facesTemp[iMax][2], facesTemp[iMax][3], 0); //pas besoin du tri ici
        faces[iMax]->construitFace(noeuds, m_numNoeuds[noeudAutre], this);
        iMax++;
      }
      else
      {
        faces[indiceFaceExiste]->ajouteElementVoisin(this);
      }
    }
    else //Faces triangles
    {
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
    } //Fin if
  } //Fin boucle faces
}

//***********************************************************************
// Nouvelle version beaucoup plus efficace avec recherche dans tableau temporaire
void ElementPrisme::construitFacesSimplifie(int &iMax, int** facesTemp, int* sommeNoeudsTemp)
{
  //3 faces a traiter de type quadrangle et 2 faces triangle
  int indiceFaceExiste(-1);
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; facesTemp[iMax][2] = m_numNoeuds[4]; facesTemp[iMax][3] = m_numNoeuds[3]; break;
      case 1: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[2]; facesTemp[iMax][2] = m_numNoeuds[5]; facesTemp[iMax][3] = m_numNoeuds[3]; break;      
      case 2: facesTemp[iMax][0] = m_numNoeuds[1]; facesTemp[iMax][1] = m_numNoeuds[2]; facesTemp[iMax][2] = m_numNoeuds[5]; facesTemp[iMax][3] = m_numNoeuds[4]; break;      
      case 3: facesTemp[iMax][0] = m_numNoeuds[0]; facesTemp[iMax][1] = m_numNoeuds[1]; facesTemp[iMax][2] = m_numNoeuds[2]; break;      
      case 4: facesTemp[iMax][0] = m_numNoeuds[3]; facesTemp[iMax][1] = m_numNoeuds[4]; facesTemp[iMax][2] = m_numNoeuds[5]; break;      
    }
    if(i<3)
    {
      sommeNoeudsTemp[iMax] = facesTemp[iMax][0] + facesTemp[iMax][1] + facesTemp[iMax][2] + facesTemp[iMax][3];
      std::sort(facesTemp[iMax], facesTemp[iMax] + 4);  //Tri des noeuds
      //Existance face ?
      indiceFaceExiste = FaceNS::rechercheFace(facesTemp[iMax], sommeNoeudsTemp[iMax], facesTemp, sommeNoeudsTemp, 4, iMax);

    }
    else
    {
      sommeNoeudsTemp[iMax] = facesTemp[iMax][0] + facesTemp[iMax][1] + facesTemp[iMax][2];
      std::sort(facesTemp[iMax],facesTemp[iMax]+3);  //Tri des noeuds
      //Existance face ?
      indiceFaceExiste = FaceNS::rechercheFace(facesTemp[iMax],sommeNoeudsTemp[iMax],facesTemp,sommeNoeudsTemp,3,iMax);
    } //Fin if
    if (indiceFaceExiste == -1)
    {
      iMax++;
    }
  }
}

//***********************************************************************

void ElementPrisme::attributFaceCommunicante(const Coord *noeuds, FaceNS **faces, const int &indiceMaxFaces, const int &nombreNoeudsInternes)
{
  int indiceFaceExiste(0);
  //Verification face 1 :
  if (m_numNoeuds[0] < nombreNoeudsInternes && m_numNoeuds[1] < nombreNoeudsInternes && m_numNoeuds[4] < nombreNoeudsInternes && m_numNoeuds[3] < nombreNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[0], m_numNoeuds[1], m_numNoeuds[4], m_numNoeuds[3]);
    if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
    {
      faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indiceFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 2 :
  if (m_numNoeuds[0] < nombreNoeudsInternes && m_numNoeuds[2] < nombreNoeudsInternes && m_numNoeuds[5] < nombreNoeudsInternes && m_numNoeuds[3] < nombreNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[0], m_numNoeuds[2], m_numNoeuds[5], m_numNoeuds[3]);
    if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
    {
      faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indiceFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 3 :
  if (m_numNoeuds[1] < nombreNoeudsInternes && m_numNoeuds[2] < nombreNoeudsInternes && m_numNoeuds[5] < nombreNoeudsInternes && m_numNoeuds[4] < nombreNoeudsInternes)
  {
    FaceQuadrangle face(m_numNoeuds[1], m_numNoeuds[2], m_numNoeuds[5], m_numNoeuds[4]);
    if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
    {
      faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indiceFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 4 :
  if (m_numNoeuds[0] < nombreNoeudsInternes && m_numNoeuds[1] < nombreNoeudsInternes && m_numNoeuds[2] < nombreNoeudsInternes)
  {
    FaceTriangle face(m_numNoeuds[0], m_numNoeuds[1], m_numNoeuds[2]);
    if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
    {
      faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indiceFaceExiste]->setEstComm(true);
    }
  }
  //Verification face 5 :
  if (m_numNoeuds[3] < nombreNoeudsInternes && m_numNoeuds[4] < nombreNoeudsInternes && m_numNoeuds[5] < nombreNoeudsInternes)
  {
    FaceTriangle face(m_numNoeuds[3], m_numNoeuds[4], m_numNoeuds[5]);
    if (face.faceExiste(faces, indiceMaxFaces, indiceFaceExiste))
    {
      faces[indiceFaceExiste]->ajouteElementVoisinLimite(this);
      faces[indiceFaceExiste]->setEstComm(true);
    }
  }
}

//***********************************************************************

int ElementPrisme::compteFaceCommunicante(vector<int*> &facesTemp, vector<int> &sommeNoeudsTemp)
{
  //3 faces a traiter de type quadrangle et 2 faces triangle
  int indiceFaceExiste(-1), nombreFacesCommunicante(0);
  int face[4], sommeNoeuds;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; face[2] = m_numNoeuds[4]; face[3] = m_numNoeuds[3]; break;
      case 1: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[2]; face[2] = m_numNoeuds[5]; face[3] = m_numNoeuds[3]; break;
      case 2: face[0] = m_numNoeuds[1]; face[1] = m_numNoeuds[2]; face[2] = m_numNoeuds[5]; face[3] = m_numNoeuds[4]; break;     
      case 3: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; face[2] = m_numNoeuds[2]; break;
      case 4: face[0] = m_numNoeuds[3]; face[1] = m_numNoeuds[4]; face[2] = m_numNoeuds[5]; break;
    }
    int iMax = sommeNoeudsTemp.size();
    if(i<3)
    {
      sommeNoeuds = face[0] + face[1] + face[2] + face[3];
      std::sort(face, face+4);
      //Recherche existance faces
      indiceFaceExiste = FaceNS::rechercheFace(face, sommeNoeuds, facesTemp, sommeNoeudsTemp, 4, iMax);
    }
    else
    {
      sommeNoeuds = face[0]+face[1]+face[2];
      std::sort(face, face+3);
      //Recherche existance faces
      indiceFaceExiste = FaceNS::rechercheFace(face,sommeNoeuds,facesTemp,sommeNoeudsTemp,3,iMax);
    }
    if (indiceFaceExiste != -1)
    {
      nombreFacesCommunicante++;
    }
  }
  return nombreFacesCommunicante;
}

//***********************************************************************
//Nouvelle version plus efficace
int ElementPrisme::compteFaceCommunicante(int &iMax, int **facesTemp, int *sommeNoeudsTemp)
{
  //3 faces a traiter de type quadrangle et 2 faces triangle
  int indiceFaceExiste(-1), nombreFacesCommunicante(0);
  int face[4], sommeNoeuds;
  for (int i = 0; i < NOMBREFACES; i++)
  {
    switch (i)
    {
      case 0: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; face[2] = m_numNoeuds[4]; face[3] = m_numNoeuds[3]; break;
      case 1: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[2]; face[2] = m_numNoeuds[5]; face[3] = m_numNoeuds[3]; break;
      case 2: face[0] = m_numNoeuds[1]; face[1] = m_numNoeuds[2]; face[2] = m_numNoeuds[5]; face[3] = m_numNoeuds[4]; break;     
      case 3: face[0] = m_numNoeuds[0]; face[1] = m_numNoeuds[1]; face[2] = m_numNoeuds[2]; break;
      case 4: face[0] = m_numNoeuds[3]; face[1] = m_numNoeuds[4]; face[2] = m_numNoeuds[5]; break;
    }
    if(i<3)
    {
      sommeNoeuds = face[0] + face[1] + face[2] + face[3];
      std::sort(face, face + 4);
      //Recherche existance faces
      indiceFaceExiste = FaceNS::rechercheFace(face, sommeNoeuds, facesTemp, sommeNoeudsTemp, 4, iMax);
    }
    else
    {
      sommeNoeuds = face[0]+face[1]+face[2];
      std::sort(face, face+3);
      //Recherche existance faces
      indiceFaceExiste = FaceNS::rechercheFace(face,sommeNoeuds,facesTemp,sommeNoeudsTemp,3,iMax);
    }
    if (indiceFaceExiste != -1)
    {
      nombreFacesCommunicante++;
    }
  }
  return nombreFacesCommunicante;
}

//***********************************************************************
