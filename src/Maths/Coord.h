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

#ifndef COORD_H
#define COORD_H

class Coord
{
public:
  Coord();
  Coord(const double &x, const double &y = 0., const double&z = 0.);
  virtual ~Coord();
  void setXYZ(const double &x, const double & y, const double &z);
  void setX(const double &x);
  void setY(const double &y);
  void setZ(const double &z);
  double getX() const;
  double getY() const;
  double getZ() const;
  double norme() const;
  double normeCarre() const;

  double scalaire(const Coord &a) const; //Produit scalaire entre vecteur et vecteur a
  Coord vectoriel(const Coord &a) const; //Produit vectoriel entre vecteur et vecteur a
  void projection(const Coord &normale, const Coord &tangente, const Coord &binormale);
  void projectionRepereAbsolu(const Coord &normale, const Coord &tangente, const Coord &binormale);

  void creeVecteur(const Coord &a, const Coord &b);
  static double produitScalaire(const Coord &v1, const Coord &v2);
  static Coord produitVectoriel(const Coord &v1, const Coord &v2);
  //void produitVectoriel(const Coord &v1, const Coord &v2);
  void changeSigne();
  void normalise();
  void afficheInfos() const;

  //Calcul le determinant de la matrice forme par les vecteurs v1,v2 et v3
  static double determinant(const Coord &v1, const Coord &v2, const Coord &v3);
  
  //Cosinus entre deux vecteurs
  static double cos(const Coord &v1, const Coord &v2);
  static Coord sin(const Coord &v1, const Coord &v2);

  //Surcharge operateurs
  Coord& operator=(const double &scalaire);
  Coord& operator*= (const double &scalaire);
  Coord& operator/= (const double &scalaire);
  Coord operator* (const double &scalaire);
  Coord operator/ (const double &scalaire);
  Coord& operator+= (const Coord &a);
  Coord& operator-= (const Coord &a);

protected:
  double m_x;
  double m_y;
  double m_z;
};

//Surcharge operateur externe a la classe car prends deux arguments
Coord operator* (const double &scalaire, const Coord &a);
Coord operator+ (const Coord &a, const Coord &b);
Coord operator- (const Coord &a, const Coord &b);

#endif // COORD_H