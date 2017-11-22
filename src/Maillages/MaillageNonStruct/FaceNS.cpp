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

#include "FaceNS.h"

using namespace std;

//***********************************************************************

FaceNS::FaceNS(){}

//***********************************************************************

FaceNS::FaceNS(const int &nombreNoeuds) :
Face(),
m_nombreNoeuds(nombreNoeuds),
m_limite(false),
m_comm(false),
m_elementGauche(0),
m_elementDroite(0)
{
  m_numNoeuds = new int[nombreNoeuds];
}

//***********************************************************************

FaceNS::~FaceNS()
{
  delete[] m_numNoeuds;
}

//***********************************************************************

ElementNS *FaceNS::getElementGauche() const
{
  return m_elementGauche;
}

//***********************************************************************

ElementNS *FaceNS::getElementDroite() const
{
  return m_elementDroite;
}

//***********************************************************************

void FaceNS::construitFace(const Coord *noeuds, const int &numNoeudAutre, ElementNS *elementVoisin)
{

  //Calcul position du centre de face
  m_position = 0.;
  for (int i = 0; i < m_nombreNoeuds; i++){ m_position += noeuds[m_numNoeuds[i]]; }
  m_position /= static_cast<double>(m_nombreNoeuds);
  //Calculs des proprietes de la face
  this->calculSurface(noeuds);
  this->calculRepere(noeuds, numNoeudAutre, elementVoisin);
}

//***********************************************************************

bool FaceNS::faceExiste(FaceNS **faces, const int &indiceMaxFaces, int &indiceFaceExiste)
{
  int faceTrouvee(1);
  //for (int i = 0; i < indiceMaxFaces; i++)
  for (int i = indiceMaxFaces-1; i>=0; i--)
  {
    //1) Sans passer par operateur de comparaison
    //-------------------------------------------
    if (faces[i]->m_sommeNumNoeuds != m_sommeNumNoeuds){continue;};
    //Verification noeud par noeud
    faceTrouvee = 1;
    for (int n = 0; n < m_nombreNoeuds; n++)
    {
      if (m_numNoeuds[n] != faces[i]->getNumNoeud(n)){ faceTrouvee = 0; break; }
    }
    if(faceTrouvee)
    {
      indiceFaceExiste = i; return true;
    }

    // //2) Avec operateur de comparaison
    // //--------------------------------
    // if (*faces[i] == *this){ indiceFaceExiste = i; return true; }
  }
  indiceFaceExiste = 0;
  return false;
}

//***********************************************************************

void FaceNS::ajouteElementVoisin(ElementNS *elementVoisin)
{
  if (m_elementGauche == 0)
  {
    m_elementGauche = elementVoisin;
  }
  else
  {
    m_elementDroite = elementVoisin;
  }
}

//***********************************************************************

void FaceNS::ajouteElementVoisinLimite(ElementNS *elementVoisin)
{
  if (m_elementGauche == 0)
  {
    //On inverse pour avoir l element de calcul a Gauche (Necessaire pour les limites)
    m_elementGauche = m_elementDroite;
    m_normale.changeSigne();
    m_binormale = Coord::produitVectoriel(m_normale, m_tangente);
  }
  m_elementDroite = elementVoisin;
  m_limite = true;
}

//***********************************************************************

void FaceNS::setEstLimite(const bool &estLimite)
{
  m_limite = estLimite;
}

//***********************************************************************

void FaceNS::setEstComm(const bool &estComm)
{
  m_comm = estComm;
}

//***********************************************************************

bool FaceNS::getEstComm() const
{
  return m_comm;
}

//***********************************************************************

int FaceNS::getNombreNoeuds() const
{
  return m_nombreNoeuds;
}

//***********************************************************************

int FaceNS::getNumNoeud(const int &numNoeud) const
{
  return m_numNoeuds[numNoeud];
}

//***********************************************************************

void FaceNS::getInfoNoeuds(int *numNoeuds, int &sommeNumNoeuds) const
{
  for(int i=0; i<m_nombreNoeuds; i++)
  {
    numNoeuds[i] = m_numNoeuds[i];
  }
  sommeNumNoeuds = m_sommeNumNoeuds;
}


//***********************************************************************

bool FaceNS::getEstLimite() const
{
  return m_limite;
}

//***********************************************************************

int FaceNS::getSommeNumNoeuds() const
{
  return m_sommeNumNoeuds;
}

//***********************************************************************

void FaceNS::afficheInfos() const
{
  cout << "----------------" << endl;
  //cout << "Face" << endl << " Noeuds : ";
  //for (int i = 0; i < m_nombreNoeuds; i++)
  //{
  //  cout << m_numNoeuds[i] << " ";
  //}
  //cout << endl;
  cout << "centre : " << m_position.getX() << " " << m_position.getY() << " " << m_position.getZ() << endl;
  //cout << " limites : " << m_limite << ", comm : " << m_comm << endl;
  //if (m_elementGauche != 0)
  //{
  //  cout << " Element gauche " << m_elementGauche->getIndice();
  //}
  //if (m_elementDroite != 0)
  //{
  //  cout << " Element Droite " << m_elementDroite->getIndice();
  //}
  //cout << endl;

  //cout << " normale : "; m_normale.afficheInfos();
  //cout << " tangente : "; m_tangente.afficheInfos();
  //cout << " binormale : "; m_binormale.afficheInfos();
}

//***********************************************************************

void FaceNS::afficheNoeuds() const
{
  cout << "Face" << " Noeuds : ";
  for (int i = 0; i < m_nombreNoeuds; i++)
  {
    cout << m_numNoeuds[i] << " ";
  }
  cout << endl;
}

//*********************************************************************

int FaceNS::rechercheFace(int *face, int &sommeNoeuds, int **tableauFaces, int *tableauSommeNoeuds, int nombreNoeuds, int &indiceMaxFaces)
{
    for (int i = indiceMaxFaces-1; i>=0; i--)
    {
      int existe = 1;
      //On utilise la somme des numero des noeuds pour faire un premier tri : permet d accelerer la recherche
      if (tableauSommeNoeuds[i] != sommeNoeuds){ continue; };
      //Verification noeud par noeud
      for (int n=0; n<nombreNoeuds; n++)
      {
        if (face[n] != tableauFaces[i][n]){ existe = 0; break; }
      }
      if (existe) return i; //Face trouvee, on renvoi son numero
    }
    return -1;  //Face non trouvee
}

//*********************************************************************

int FaceNS::rechercheFace(int *face, int &sommeNoeuds, vector<int*> tableauFaces, vector<int> tableauSommeNoeuds, int nombreNoeuds, int &indiceMaxFaces)
{
    for (int i = indiceMaxFaces-1; i>=0; i--)
    {
      int existe = 1;
      //On utilise la somme des numero des noeuds pour faire un premier tri : permet d accelerer la recherche
      if (tableauSommeNoeuds[i] != sommeNoeuds){ continue; };
      //Verification noeud par noeud
      for (int n=0; n<nombreNoeuds; n++)
      {
        if (face[n] != tableauFaces[i][n]){ existe = 0; break; }
      }
      if (existe) return i; //Face trouvee, on renvoi son numero
    }
        // cout << "AHAHAHAHAHA nouveau" << endl;
    return -1;  //Face non trouvee
}

//*********************************************************************
//Surcharge operateur externe a la classe car prends deux arguments

bool operator==(const FaceNS &a, const FaceNS &b)
{
  //Sortie directe si m_sommeNumNoeuds differents
  if (a.getSommeNumNoeuds() != b.getSommeNumNoeuds()){ return false; }
  //Sortie directe si nombre noeuds differents
  if (a.getNombreNoeuds() != b.getNombreNoeuds()){ return false; }
  //Verification noeud par noeud
  for (int i = 0; i < a.getNombreNoeuds(); i++)
  {
    if (a.getNumNoeud(i) != b.getNumNoeud(i)){ return false; }
  }
  //Par defaut : il sont identiques
  return true;
}

//*********************************************************************

bool operator!=(const FaceNS &a, const FaceNS &b)
{
  return !(a == b);
}

//*********************************************************************