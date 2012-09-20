#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#define MINI(X,Y)         ( (X) < (Y) ? (X) : (Y) )
#define MAXI(X,Y)         ( (X) > (Y) ? (X) : (Y) )

// To switch on the support for complex-valued matrices, you need to define the macro LA_COMPLEX_SUPPORT in your application before including the Lapack++ header files.
#define LA_COMPLEX_SUPPORT

#include <algorithm>
#include <complex>
#include <ctime>
#include <fstream>
#include <gmd.h>
#include <iostream>
#include <sstream>
#include <vector>

// LAPACKPP
#include "lafnames.h"       /* macros for LAPACK++ filenames */
#include LA_GEN_MAT_COMPLEX_H
#include LA_VECTOR_DOUBLE_H
#include "blaspp.h"
#include LA_SOLVE_DOUBLE_H
#include LA_GENERATE_MAT_DOUBLE_H
#include LA_EXCEPTION_H
#include LA_UTIL_H
#include LA_PREFS_H 
#include "lasvd.h"


#include "Particle.hpp"

template <typename T>
bool isIn(T x, std::vector<T> v)
{
  unsigned int i=0;
  while(i<v.size() && x != v[i])
    ++i;
  return (i != v.size());
}


LaVectorComplex genMatComplexToVectorComplex(LaGenMatComplex & m)
{
  LaVectorComplex res(m.rows());
  for(unsigned int i=0;i<m.rows();++i)
    res(i)=m(i,0);
  return res;
}



std::vector<std::vector<Particle> > mesh(unsigned int x, unsigned int y)
{
  std::cerr << "Creating the mesh..." << std::endl;
  std::vector<std::vector<Particle> > res;
  for(unsigned int i=0;i<x-1;++i)
  {
//    fprintf(stderr,"%lf%%\r",(double)(i*100)/(x*y));
    std::vector<Particle> row;
    for(unsigned int j=0;j<y-1;++j)
    {
      //res.push_back(Particle( ((double)i +1)/((double)x) ,  ((double)j+1)/((double)y)  ) );
      row.push_back(Particle( ((double)i +1)/((double)x) ,  ((double)j+1)/((double)y)  ) );
    }
    res.push_back(row);
  }
  std::cerr << " done." << std::endl;
  return res;
}

/*
// Mesh of all sequences of a given set of length N
std::vector<std::vector<int> > mesh3Q(int N)
{
  std::vector<std::vector<int> > res;
  
  std::vector<int> seq;
  int index=0; 
  for(int i=0;i<N;++i)
  {
    seq.push_back(0);
  }
  res.push_back(seq);
  bool back=false;
  while(index<N)
  {
    int j=0;
    while(j<=index)
    {
      if(seq[j]==2)
      {
        while(seq[j]==2 && j<=index)
        {
          ++j;
        }

      }
      else
      {
        ++seq[j];
        for(int k=0;k<j;++k)
          seq[k]=0;
        res.push_back(seq);
        j=0;
      }
    }
    ++index;
    if(index<N)
    {
      ++seq[index];
      for(int k=0;k<index;++k)
        seq[k]=0;
      res.push_back(seq);
    }
  }
  
  return res;
}



*/

/*                                                             */
/* Given a int * representing the periodic sequence */
/* in base 3, return the decimal value                         */
/*                                                             */
double decValue(int * v, int n)
{
  double res = 0;
  for(int i=0;i<n;++i)
    res+=v[i]*pow(3,i);
  res = res/(pow(3,n)-1);
  return res;
}




/*                                    */
/* Creates a mesh with a step of 3^-N */
/*                                    */

std::vector<Particle> meshPower3(int N)
{
  std::vector<Particle> res;
  double step = pow(3,-N);
  
//  std::cout << "N " << N << std::endl;
//  std::cout << step << std::endl;
  
  double q = step;
  double p = 0.5;
  
  while(q<1)
  {
//    std::cout << "(" << q << "," << p << ")" << std::endl;
    res.push_back(Particle(q,p));
    q += step;
  }
  
  q=step;
  p+=step;
  while(q<1)
  {
    
//    std::cout << "(" << q << "," << p << ")" << std::endl;
    res.push_back(Particle(q,p));
    q += step;
  }
  
  return res;
}


/*                                             */
/*  Calculate the adjoint of a LaGenMatComplex */
/*                                             */

LaGenMatComplex adjoint(LaGenMatComplex & m)
{
  LaGenMatComplex res(m.rows(),m.cols());
  
  for(int i=0;i<m.rows();++i)
    for(int j=0;j<m.cols();++j)
      {
        std::complex<double> z(m(i,j).r,m(i,j).i);
        LaComplex zz(std::conj(z));
        res(j,i)=zz;
      }

  return res;
}


/*                                    */
/* Constructs the FN Matrix of size N */
/*                                    */

LaGenMatComplex FNMatrix(unsigned int N)
{
  LaGenMatComplex res(N,N);
  double alpha = 1/sqrt(N);
  double PI = 3.141592653589793238462643383279;
  
  LaComplex I(0,1);
  
  for(int m=0;m<N;++m)
  {
    fprintf(stderr,"%lf%%\r",(double)(m*100)/N);
    for(int n=0;n<N;++n)
    {
//      std::complex<double> z(0,-2*PI/N*(m+0.5)*(n+0.5));
      std::complex<double> z(0,-2*PI/N*(m+0.5)*(n+0.5));
      LaComplex zz(alpha*exp(z));
      res(m,n) = zz;
    }
  }

  return res;
}




/*                                        */
/*  Construct the diagonal block matrix : */
/*                                        */
/*  FN/3           0          0           */
/*   0     sqrt(rho)*FN/3     0           */
/*   0             0         FN/3         */
/*                                        */
LaGenMatComplex diagBlock(LaGenMatComplex & FN3, double rho)
{
  int m = FN3.rows();
  int n = FN3.cols();
  LaGenMatComplex res(LaGenMatComplex::zeros(3*m,3*n));
  
  if(0 <= rho && rho <=1)
  {
    double sqrtRho = sqrt(rho);
    
    for(int i=0;i<m;++i)
      for(int j=0;j<n;++j)
        res(i,j) = FN3(i,j);
    
    for(int i=0;i<m;++i)
      for(int j=0;j<n;++j)
        res(m+i,n+j) = LaComplex(FN3(i,j).r*sqrtRho, FN3(i,j).i*sqrtRho);
    
    for(int i=0;i<m;++i)
      for(int j=0;j<n;++j)
        res(m*2+i,n*2+j) = FN3(i,j);
  }
  else
    std::cerr << "Error in diagBlock(LaGenMatComplex & FN3, double rho) : rho is not in [0,1]. " << std::endl;
 
  return res;
}

LaComplex conjugate(LaComplex z)
{
  return LaComplex(z.real(),-z.imag());
}


LaComplex conjugate(std::complex<double> & z)
{
  return LaComplex(z.real(),-z.imag());
}

/*                                       */
/* Complex scalar product                */
/* vectors are here LaGenMatComplex(N,1) */
/*                                       */

LaComplex complexScalarProd(LaGenMatComplex a, LaGenMatComplex b)
{
  LaComplex res(0,0);
  double ar,ai,br,bi;
  double re,im;
  for(unsigned int i=0;i<a.rows();++i)
  {
/*
    ar = a(i,0).r;
    ai = a(i,0).i;
    br = b(i,0).r;
    bi = b(i,0).i;
    re = ar*br + ai*bi;
    im = ar*bi - br*ai;
*/    
    LaComplex aa(a(i,0));
    LaComplex bb(b(i,0));

//    re = ar*br - ai*bi;
//    im = ar*bi + br*ai;
    res = res + conjugate(aa)*bb;
//    res = LaComplex(res.real() + re, res.imag() + im);

//    conjugate(a(i,1)
  }
  return res;
}

/*                                                     */
/* Return the components of a vector in the eigenbasis */
/*                                                     */

LaGenMatComplex toEigenBasis(LaGenMatComplex & psi, LaGenMatComplex & eigVecL, LaGenMatComplex & eigVecR)
{
  LaGenMatComplex res(eigVecL.rows(),1);
  LaGenMatComplex vR(eigVecR.rows(),1);
  LaGenMatComplex vL(eigVecL.rows(),1);

  double epsilon = 0.00000001;

  for(unsigned int i=0;i<eigVecL.cols();++i)
  {
    for(unsigned int k=0;k<eigVecR.cols();++k)
      vR(k,0) = eigVecR(k,i);

    for(unsigned int k=0;k<eigVecL.cols();++k)
      vL(k,0) = eigVecL(k,i);
  
    if(norm(complexScalarProd(vL,vR)) > epsilon)
      res(i,0) = (complexScalarProd(vL,psi)/complexScalarProd(vL,vR));
    else
      std::cerr << "ERROR in toEigenBasis : " << i << " : " << complexScalarProd(vL,vR) << std::endl;
  }
  return res;
}


/*                                         */
/* Sort left eigenvalues in the same order */
/* than the right eigenvectors             */
/*                                         */

LaGenMatComplex sortLeftEig(LaGenMatComplex & eigVecR, LaGenMatComplex & eigVecL)
{
  LaGenMatComplex res(eigVecR.rows(),eigVecR.cols());
  LaGenMatComplex vR = eigVecR.col(0);
  LaGenMatComplex vL = eigVecL.col(0);
  unsigned int j=0;
  std::vector<unsigned int> indexes;
  bool cont = true;
  double epsilon = 0.000000001;
  unsigned int N = eigVecR.cols();

//  std::cerr << "kokoDAA" << std::endl;
  for(unsigned int i=0;i<N;++i)
  {
    fprintf(stderr,"%lf%%\r",(double)(i*100)/N);
    for(unsigned int k=0;k<eigVecR.cols();++k)
      vR(k,0) = eigVecR(k,i);

//    vR = eigVecR.col(i).copy();
//  std::cerr << "kokoDAA1" << std::endl;
   for(unsigned int k=0;k<eigVecR.cols();++k)
      vL(k,0) = eigVecL(k,0);

//    vL = eigVecL.col(0).copy();
//  std::cerr << "kokoDAA2" << std::endl;
    j=0;
    cont = true;
    while(j<eigVecL.cols() && cont)
    {
      for(unsigned int k=0;k<eigVecR.cols();++k)
        vL(k,0) = eigVecL(k,j);

//      vL = eigVecL.col(j).copy();
      if(norm(complexScalarProd(vR,vL)) > epsilon)
      {
        if(!isIn(j,indexes))
        {
          cont = false;
          for(unsigned int k=0;k<eigVecL.cols();++k)
          {
//            std::cerr << "kokoiFOR" << std::endl;
            res(k,i) = eigVecL(k,j);
          }
        }
        else
        {
          std::cerr << "Probleme : i:" << i << " , j:" << j << " , ps:" << complexScalarProd(vR,vL) << std::endl;
          exit(0);
        }
      }
      ++j;
    }
  }
  return res;
}

/*                           */
/* Power each component to n */
/*                           */
LaVectorDouble powerVector(LaVectorDouble & v, unsigned int n )
{
  LaVectorDouble res(v.size());
  for(unsigned int i=0;i<v.size();++i)
    res(i) = pow(v(i),n);
  return res;
}

/*                                                            */
/* Return the expecting value of a vector for the position    */
/* operator, after n iterations :  <x> = <v Un | x | Un v>    */
/* v must be given in the eigenbasis representation           */
/* The order of the eigenvalues determines the representation */
/*                                                            */
LaComplex expectingValue(LaGenMatComplex & v, LaGenMatComplex & eigVecR, LaGenMatComplex & eigValR, unsigned int n)
{
  std::cout << "Calculating the expecting value..." << std::endl;
  LaComplex res(0,0);
  LaComplex cili,cjlj,psiik,psijk;
  unsigned int N=v.rows();
  double qk=0;
  LaVectorComplex vec(N);
  
  LaComplex vi;
  // calculate the c_m*(lambda_m)^n
  for(unsigned int i=0;i<N;++i)
  {
    vi = v(i,0);
//    std::cout << vi << std::endl;
    vec(i) = LaComplex(pow(std::complex<double>(eigValR(i,0).r,eigValR(i,0).i),n))*vi;
//    if(vec(i).r > 0.000001)
//      std::cout << vec(i) << std::endl;
  }
  LaComplex restmp(0,0); 
  for(unsigned int k=0;k<N;++k)
  {
    fprintf(stderr,"%lf%%\r",(double)(k*100)/N);
    qk = (k+0.5)/N;
    for(unsigned int i=0;i<N;++i)
    {
      cili  = conjugate(vec(i));
      psiik = conjugate(eigVecR(k,i));
      //std::cout << "cili : " << cili << "   psiik : " << psiik << std::endl;
      for(unsigned int j=0;j<N;++j)
      {
        cjlj = vec(j);
        psijk = eigVecR(k,j);
//        std::cout << "cili*cjlj*psiik*psijk : " << cili*cjlj*psiik*psijk << std::endl;
//        restmp = 
        res +=  cili*cjlj*psiik*psijk;
      }
      std::cout << "res : " << res << std::endl << std::endl;
    }
    res *= qk;
  }
/*    



  {
    lici = conjugate(vec(i));
//    ci = conjugate(v(i,0));
//    for(unsigned int j=0;j<v.rows();++j)
    {
    int j=i;
//      ci = conjugate(v(i,0));
      cj = v(j,0);
//      li = conjugate(pow(std::complex<double>(eigValR(i,0).r,eigValR(i,0).i),n));
//      lj = pow(std::complex<double>(eigValR(j,0).r,eigValR(j,0).i),n);
      lj = eigV(j);
//      scal = complexScalarProd(eigVecR.col(i),eigVecR.col(j));i
//      scal = eigVecScal(i,j);

      std::cout << i << "/" << j << std::endl;
      std::cout << "   ci   : " << ci << std::endl;
      std::cout << "   cj   : " << cj << std::endl;
      std::cout << "   li   : " << li << std::endl;
      std::cout << "   lj   : " << lj << std::endl;
      std::cout << "   scal : " << scal << std::endl;



      res = res + ci*cj*li*lj;
//      res+= conjugate(v(i,0)) * v(j,0) * conjugate(pow(eigValR(i,0),n)) * pow(eigValR(j),n) * complexScalarProduct(eigVecR.col(i),eigVecR.col(j));
    }
  }
*/
  return res;
}


/*                                                            */
/* Returns the standard deviation of the vector v (which is   */
/* actually a matrix with 1 column) for the position operator */
/* deltaX = sqrt( <x^2> - <x>^2 )                             */
/*        = sqrt( <v|x^2|v> - <v|x|v>^2)                      */
/* v must be given in the eigenbasis representation           */
/* The order of the eigenvalues determines the representation */
/*                                                            */

double incertitudeX(std::vector<Particle> & mesh, LaVectorComplex & eigValSort, LaGenMatComplex & v)
{
  double res=0;
  for(unsigned int i=0;i<v.rows();++i)
  {
  
  }
  return res;
}




/*                                                           */
/* Returns average of the husimi representation on n vectors */
/* (n is the number of columns of the matrix v)              */
/*                                                           */

LaVectorDouble husimi(std::vector<std::vector<Particle> > & mesh, LaGenMatComplex & v )
{
  std::cout << "Calculating the average Husimi representation..." << std::endl;
  
  LaVectorDouble hus(mesh.size()*mesh[0].size());
  int N=v.rows();
  int n=v.cols(); // (the n sorted eigenvectors)
  
  double PI = 3.141592653589793238462643383279;
  LaComplex I(0,1);
  double q,p,sum;
  LaComplex sum1(0,0);
  LaComplex sum2(0,0);
  double ekm;
  std::complex<double> z;
  LaComplex zz;
  LaComplex vkr;
 
 
//std::cout << "mesh.size()" << mesh.size() << "   " << mesh[0].size() << std::endl;
//std::cout << "matrix size()" << v.size(0) << "   " << v.size(1) << std::endl;

  for(unsigned int partx=0;partx<mesh.size();++partx)
  {
    fprintf(stderr,"%lf%%\r",(double)(partx*100)/mesh.size());
    for(unsigned int party=0;party<mesh[0].size();++party)
    {
      q=mesh[partx][party].getQ();
      p=mesh[partx][party].getP();
//*
      sum=0;
      for(int r=0;r<n;++r)
      {
        sum1 = LaComplex(0,0);
        for(int k=0;k<N;++k)
        {
          sum2 = LaComplex(0,0);
          for(int m=-3;m<=3;++m)
          {
            ekm=(k-0.5)/N + m;
            z = std::complex<double>(-PI*N*pow(ekm - q,2) , 2*PI*N*(ekm - q/2)*p - PI*m);
            zz = exp(z);
            sum2 += zz;
          }
          vkr = v(k,r);
          sum1 += vkr*sum2;
        }
        sum += pow(abs(sum1),2);
      }
      sum = sum/n;
//*/
/*
      if(10<partx && partx<15 && 30<party && party<35)
        hus(partx*mesh[0].size()+party) = 1;
      else
        hus(partx*mesh[0].size()+party) = 0;
//*/
//      hus(partx*mesh[0].size()+party) = p;
    }
  }
  return hus;
}






/*                                                         */
/*  Returns the average of the amplitude of the eigensates */
/*                                                         */

std::vector<double> amplitude(LaGenMatComplex v )
{
  std::cout << "Calculating the amplitude..." << std::endl;
  
  std::vector<double> ampl;
  int N=v.rows();
  int n=v.cols(); // (the n sorted eigenvectors)
  double sum;
  for(unsigned int i=0;i<N;++i)
  {
    sum = 0;
    for(unsigned int j=0;j<n;++j)
    {
      sum += v(i,j).r*v(i,j).r + v(i,j).i*v(i,j).i;
    }
    sum = sum/n;
  }
  
  return ampl;
}





/*                                                    */
/* Read and create a LaVectorComplex with 2 files :   */
/* The first one contains the real part of the vector */
/* The second one contains the imaginary part         */
/*                                                    */
LaVectorComplex readReImVector(char * fileNameRe, char * fileNameIm)
{
  std::cerr << "Reading and creating vector..." << std::endl;

  // Real part
  std::vector<double> re;
  std::fstream fileRe(fileNameRe,std::ios::in);
  if(fileRe)
  {
    std::string line;
    double x;
    std::stringstream sline;
    
    while(getline(fileRe,line))
    {
      sline.clear();
      sline.str("");
      sline << line;
      //while(sline >> x)
      sline >> x;
        re.push_back(x);
    }
    fileRe.close();
  }
  else
    std::cerr << "Error in readReImVector(char * fileNameRe, char * fileNameIm) : fileRe is not here" << std::endl;
  
  
  // Imaginary part
  std::vector<double> im;
  std::fstream fileIm(fileNameIm,std::ios::in);
  if(fileIm)
  {
    std::string line;
    double x;
    std::stringstream sline;
    
    while(getline(fileIm,line))
    {
      sline.clear();
      sline.str("");
      sline << line;
      //while(sline >> x)
      sline >> x;
        im.push_back(x);
    }
    fileIm.close();
  }
  else
    std::cerr << "Error in readReImVector(char * fileNameRe, char * fileNameIm) : fileIm is not here" << std::endl;
  
  LaVectorComplex res(re.size());
  
  for(unsigned int i=0;i<re.size();++i)
  {
    LaComplex z(re[i],im[i]);
    res(i) = z;
  }
  
  std::cerr << " done." << std::endl;
  return res;
}


/*                                                    */
/* Read and create a LaGenMatComplex with 2 files :   */
/* The first one contains the real part of the matrix */
/* The second one contains the imaginary part         */
/*                                                    */
LaGenMatComplex readReImMatrix(char * fileNameRe, char * fileNameIm)
{
  std::cerr << "Reading and creating matrix..." << std::endl;

  // Real part
  std::vector< std::vector<double> > matRe;
  std::fstream fileRe(fileNameRe,std::ios::in);
  if(fileRe)
  {

    std::string line;
    double x;
    std::stringstream sline;
    
    while(getline(fileRe,line))
    {
      std::vector<double> row;
      sline.clear();
      sline.str("");
      sline << line;
      while(sline >> x)
        row.push_back(x);
      matRe.push_back(row);
    }
    fileRe.close();
  }
  else
    std::cerr << "Error in readReImVector(char * fileNameRe, char * fileNameIm) : fileRe is not here" << std::endl;
  LaGenMatDouble re(matRe.size(),matRe[0].size());
  for(unsigned int i=0;i<matRe.size();++i)
    for(unsigned int j=0;j<matRe[0].size();++j)
      re(i,j)=matRe[i][j];

  
  // Imaginary part
  std::vector< std::vector<double> > matIm;
  std::fstream fileIm(fileNameIm,std::ios::in);
  if(fileIm)
  {
    std::string line;
    double x;
    std::stringstream sline;
    
    while(getline(fileIm,line))
    {
      std::vector<double> row;
      sline.clear();
      sline.str("");
      sline << line;
      while(sline >> x)
        row.push_back(x);
      matIm.push_back(row);
    }
    fileIm.close();
  }
  else
    std::cerr << "Error in readReImVector(char * fileNameRe, char * fileNameIm) : fileIm is not here" << std::endl;
  LaGenMatDouble im(matIm.size(),matIm[0].size());
  for(unsigned int i=0;i<matIm.size();++i)
    for(unsigned int j=0;j<matIm[0].size();++j)
      im(i,j)=matIm[i][j];

  LaGenMatComplex res(re.rows(),re.cols());
  for(unsigned int i=0;i<res.rows();++i)
    for(unsigned int j=0;j<res.cols();++j)
    {
      res(i,j).r=re(i,j);
      res(i,j).i=im(i,j);
    }
  std::cerr << " done." << std::endl;
  return res;


}

void writeLaVectorDouble(LaVectorDouble & v, char * fileName)
{
  std::cerr << "Storing..." << std::endl;
  std::fstream fileOut(fileName,std::ios::out);
  if(fileOut)
  {
    for(unsigned int i=0;i<v.size();++i)
      fileOut << v(i) << std::endl;
    fileOut.close();
  }
  else
    std::cerr << "writeLaVectorDouble(LaVectorDouble & v, char * fileName) : can't create fileName." << std::endl;
}



void writeVector(std::vector<unsigned int> &v, char * fileName)
{
  std::cerr << "Storing..." << std::endl;
  std::fstream fileOut(fileName,std::ios::out);
  if(fileOut)
  {
    for(unsigned int i=0;i<v.size();++i)
      fileOut << v[i] << std::endl;
    fileOut.close();
  }
  else
    std::cerr << "writeVector(std::vector<unsigned int> & v, char * fileName) : can't create fileName." << std::endl;
}

std::vector<unsigned int> readVector(char * fileName)
{
  std::vector<unsigned int> res;

  std::fstream fileIn(fileName,std::ios::in);
  if(fileIn)
  {
    std::string line;
    unsigned int x;
    std::stringstream sline;

    while(getline(fileIn,line))
    {
      sline.clear();
      sline.str("");
      sline << line;
      sline >> x;
      res.push_back(x);
    }
   fileIn.close();

  }
  else
    std::cerr << "readVector(char * fileName) : can't create fileName." << std::endl;

  return res;
}


/*                                                    */
/* Write a LaVectorComplex with 2 files :             */
/* The first one contains the real part of the vector */
/* The second one contains the imaginary part         */
/*                                                    */
void writeReImVector(LaVectorComplex & V, char * fileNameRe, char * fileNameIm)
{
  std::cerr << "Storing..." << std::endl;

  // Real part
  std::fstream fileOut1(fileNameRe,std::ios::out);
  if(fileOut1)
  {
    for(unsigned int i=0;i<V.size();++i)
      fileOut1 << V(i).r << std::endl;
    fileOut1.close();
  }
  else
    std::cerr << "writeReImVector(LaVectorComplex & V, char * fileNameRe, char * fileNameIm) : can't create fileNameRe." << std::endl;

  // Imaginary part
  std::fstream fileOut2(fileNameIm,std::ios::out);
  if(fileOut2)
  {
    fileOut2.precision(15);
    for(unsigned int i=0;i<V.size();++i)
      fileOut2 << V(i).i << std::endl;
    fileOut2.close();
  }
  else
    std::cerr << "writeReImVector(LaVectorComplex & V, char * fileNameRe, char * fileNameIm) : can't create fileNameIm." << std::endl;
}

/*                                                    */
/* Write a LaGenMatComplex in 2 files :               */
/* The first one contains the real part of the matrix */
/* The second one contains the imaginary part         */
/*                                                    */
void writeReImMatrix(LaGenMatComplex & M, char * fileNameRe, char * fileNameIm)
{
  std::cerr << "Storing..." << std::endl;
  
  // Real part
  std::fstream fileOut1(fileNameRe,std::ios::out);
  if(fileOut1)
  {
    fileOut1.precision(16);
    for(unsigned int i=0;i<M.rows();++i)
    {
      for(unsigned int j=0;j<M.cols();++j)
        fileOut1 << M(i,j).r << " ";
      fileOut1 << std::endl;
    }
    fileOut1.close();
  }
  else
    std::cerr << "Error in writeReImMatrix(LaGenMatComplex & M, char * fileNameRe, char * fileNameIm) : can't create fileNameRe." << std::endl;
  
  // Imaginary part
  std::fstream fileOut2(fileNameIm,std::ios::out);
  if(fileOut2)
  {
    fileOut2.precision(16);
    for(unsigned int i=0;i<M.rows();++i)
    {
      for(unsigned int j=0;j<M.cols();++j)
        fileOut2 << M(i,j).i << " ";
      fileOut2 << std::endl;
    }
    fileOut2.close();
  }
  else
    std::cerr << "Error in writeReImMatrix(LaGenMatComplex & M, char * fileNameRe, char * fileNameIm) : can't create fileNameIm." << std::endl;

}



/*                                                */
/* Given a matrix whose columns are eigen vectors */
/* of size M and a vector representing the        */
/* corresponding eigen values, return a           */
/* LaGenMatComplex of size MxN, composed of the   */
/* N eigen vectors corresponding to the Nth       */
/* highest moduli eigen values.                   */
/*                                                */
//LaGenMatComplex sortHighest(LaGenMatComplex & eigVect, LaVectorComplex & eigVal, unsigned int N)
std::pair<LaGenMatComplex,LaVectorComplex> sortHighest(LaGenMatComplex & eigVect, LaVectorComplex & eigVal, unsigned int N)
{
  LaGenMatComplex res(eigVect.rows(),N);
  LaVectorComplex res2(N);
  
  if(N <= eigVect.cols())
  {
    std::vector<unsigned int> index;
    double max=-1;
    int maxIndex=-1;
    unsigned int i=0;
    double sqNorm=0;
    for(unsigned int t=0;t<N;++t)    
    {
      maxIndex = -1;
      max = -1;
      i = 0;
      while(i<eigVal.size())
      {
        if(!isIn(i,index))
        {
          sqNorm = eigVal(i).r*eigVal(i).r + eigVal(i).i*eigVal(i).i;
          if(max < sqNorm)
          {
            max = sqNorm;
            maxIndex = i;
          }
        }
        ++i;
      }
      index.push_back(maxIndex);
      
      for(unsigned int j=0;j<eigVect.rows();++j)
        res(j,t) = eigVect(j,maxIndex);
      res2(t) = eigVal(maxIndex);
    }
  }
  else
    std::cerr << "Error in sortHighest(LaGenMatComplex eigVect, LaVectorComplex eigVal, unsigned int N) : N > nb columns (" << N << ">" << eigVect.cols() << ")" << std::endl;
  return std::pair<LaGenMatComplex,LaVectorComplex>(res,res2);
}



/*                                                */
/* Given a matrix whose columns are eigen vectors */
/* of size M and a vector representing the        */
/* corresponding eigen values, return a           */
/* LaGenMatComplex of size MxN, composed of the   */
/* N eigen vectors corresponding to the Nth       */
/* lowest moduli eigen values.                    */
/*                                                */
//LaGenMatComplex sortLowest(LaGenMatComplex & eigVect, LaVectorComplex & eigVal, unsigned int N)
std::pair<LaGenMatComplex,LaVectorComplex> sortLowest(LaGenMatComplex & eigVect, LaVectorComplex & eigVal, unsigned int N)
{
    LaGenMatComplex res(eigVect.rows(),N);
    LaVectorComplex res2(N);
  
  
//  LaGenMatComplex res(eigVect.rows(),N);
  if(N <= eigVect.cols())
  {
    std::vector<unsigned int> index;
    double min=1.1;
    int minIndex=-1;
    unsigned int i=0;
    double sqNorm=0;
    for(unsigned int t=0;t<N;++t)    
    {
      minIndex = -1;
      min = 1.1;
      i = 0;
      while(i<eigVal.size())
      {
        if(!isIn(i,index))
        {
          sqNorm = eigVal(i).r*eigVal(i).r + eigVal(i).i*eigVal(i).i;
          if(min > sqNorm)
          {
            min = sqNorm;
            minIndex = i;
          }
        }
        ++i;
      }
      index.push_back(minIndex);
      
      for(unsigned int j=0;j<eigVect.rows();++j)
        res(j,t) = eigVect(j,minIndex);
      res2(t) = eigVal(minIndex);
    }
  }
  else
    std::cerr << "Error in sortLowest(LaGenMatComplex eigVect, LaVectorComplex eigVal, unsigned int N) : N > nb columns (" << N << ">" << eigVect.cols() << ")" << std::endl;
  return std::pair<LaGenMatComplex,LaVectorComplex>(res,res2);
}


/*                                             */
/* Given a vector of eigenvalues, return the N */
/* highest moduli                              */
/*                                             */

LaVectorDouble sortHighest(LaVectorComplex & eigVal, unsigned int N)
{
  std::vector<double> v;
  for(unsigned int i=0;i<eigVal.size();++i)
    v.push_back(sqrt(eigVal(i).r*eigVal(i).r + eigVal(i).i*eigVal(i).i));

  std::sort(v.begin(),v.end());
  std::reverse(v.begin(),v.end());
  
  LaVectorDouble res(N);
  for(unsigned int i=0;i<N;++i)
    res(i) = v[i];

  return res;
}


/*                                             */
/* Given a vector of eigenvalues, return the N */
/* lowest moduli                               */
/*                                             */

LaVectorDouble sortLowest(LaVectorComplex & eigVal, unsigned int N)
{
  std::vector<double> v;
  for(unsigned int i=0;i<eigVal.size();++i)
    v.push_back(sqrt(eigVal(i).r*eigVal(i).r + eigVal(i).i*eigVal(i).i));

  std::sort(v.begin(),v.end());
  
  LaVectorDouble res(N);
  for(unsigned int i=0;i<N;++i)
    res(i) = v[i];

  return res;
}


/*                                                  */
/* Given a mesh and a vector of eigenvalues, return */
/* the histogram of the spectrum of the eigenvalue  */
/* moduli : for an x, it calculates the number of   */
/* eigenvalues whose moduli is < x                  */
/*                                                  */

LaVectorDouble spectrumModuliCumulated(std::vector<std::vector<Particle> > & mesh, LaVectorComplex & eigVal)
{
  LaVectorDouble res(mesh[0].size()+1);
  
  LaVectorDouble sortedModuli = sortLowest(eigVal,eigVal.size());

  unsigned int k=0;
  unsigned int nb = 0;

  unsigned int i=0;
  unsigned int j=0;

  for(unsigned int i=0;i<mesh.size();++i)
  {
    while(k<sortedModuli.size() && sortedModuli(k)<mesh[i][0].getQ())
    {
      ++nb;
      ++k;
    }
    res(i) = nb;
  }
  return res;
}


/*                                                  */
/* Return the spectrum of the eigenvalue moduli     */
/* The result is sorted from the highest moduli to  */
/* the lowest                                       */
/*                                                  */
LaVectorDouble spectrumModuli(LaVectorComplex & eigVal)
{
  return sortHighest(eigVal,eigVal.size());
}


/*                                                             */
/* Return the number of eigenvalues whose moduli is equal to x */
/*                                                             */
unsigned int spectrumModuliEqual(double x, LaVectorComplex & eigVal)
{
  unsigned int res=0;
  for(unsigned int i=0;i<eigVal.size();++i)
    if( sqrt(eigVal(i).r*eigVal(i).r + eigVal(i).i*eigVal(i).i) == x )
      ++res;
  return res;
}



/*                                                             */
/* Given the matrix whose columns are eigenvalues of Utilda,   */
/* a vector containing the corresponding eigenvalues,          */
/* the value around which we want to find the eigenvalues, and */
/* the number of eigenvectors we want to have,                 */
/* returns the average of the amplitude for these eigenvectors */
/*                                                             */

LaVectorDouble averageAmplitude(LaGenMatComplex & eigVec, LaVectorComplex & eigVal, double value, unsigned int n)
{
  std::cerr << "Calculating the average of amplitude..." << std::endl;
  LaVectorDouble res(eigVec.rows());

  std::cerr << "Sorting the vectors..." << std::endl;
  std::pair<LaGenMatComplex,LaVectorComplex> S = sortLowest(eigVec,eigVal,eigVec.rows());
  
  // Check if the parameters are coherent
  if(value >= 0 && value <= 1 && n <= eigVec.cols())
  {

    // Return the index of the eigenvalue whose moduli is immediatly greater than value
    unsigned int index=0;
    while(index<S.second.size() && sqrt(S.second(index).r*S.second(index).r + S.second(index).i*S.second(index).i) < value)
      ++index;


    unsigned int indexBegin, indexEnd;
    indexBegin = MAXI(0,index-n/2);
    indexEnd = MINI(eigVal.size(),index+n/2);
    
    if(indexEnd - indexBegin != n)
      std::cerr << "Warning : the average is done on " << indexBegin - indexEnd << " values, instead of " << n << "." << std::endl;
    
    // Average
    double sum=0;
    for(unsigned int i=0;i<eigVec.rows();++i)
    {
      sum=0;
      for(unsigned int j=indexBegin;j<indexEnd;++j)
      {
        sum += S.first(i,j).r*S.first(i,j).r + S.first(i,j).i*S.first(i,j).i;
//        std::cout << S.first(i,j) << " ";
      }
      std::cout << std::endl << sum << std::endl;
      res(i) = sum/(double)n;
    }

    std::cerr << " done." << std::endl;
/*
// for testing, adel    
    std::cerr << "[" << value << "] -> index : " << index << std::endl;
    
    if(value == 0)
    {
      for(unsigned int i=0;i<S.second.size();++i)
        std::cout << i << " " << sqrt(S.second(i).r*S.second(i).r + S.second(i).i*S.second(i).i) << std::endl;
    }
    
    
    index = 0;
    while(index<S.second.size() && sqrt(S.second(index).r*S.second(index).r + S.second(index).i*S.second(index).i) < 0.4)
      ++index;
    
    std::cerr << "[" << 0.4 << "] -> index : " << index << std::endl;
    
    index = 0;
    while(index<S.second.size() && sqrt(S.second(index).r*S.second(index).r + S.second(index).i*S.second(index).i) < 0.5)
      ++index;
    
    std::cerr << "[" << 0.5 << "] -> index : " << index << std::endl;
    
    index = 0;
    while(index<S.second.size() && sqrt(S.second(index).r*S.second(index).r + S.second(index).i*S.second(index).i) < 0.8)
      ++index;
    
    std::cerr << "[" << 0.8 << "] -> index : " << index << std::endl;
    
    
    ;
*/
  }
  else
    std::cerr << "Error in averageAmplitude(LaGenMatComplex & eigVec, LaVectorComplex & eigVal, double value, unsigned int n) : Check the parameters." << std::endl;

  return res;
}


/*                                   */
/* Return the coherent state : |q,p> */
/*                                   */

LaVectorComplex coherentState(double q, double p, unsigned int N)
{
  LaVectorComplex res(N);
  double PI=3.14159265358979323;
  LaComplex sum2(0,0);
  std::complex<double> z(0,0);
  LaComplex zz(0,0);
  double ekm=0;
  
  
  for(unsigned int k=0;k<N;++k)
  {
    sum2 = LaComplex(0,0);
    for(int m=-3;m<=3;++m)
    {
      ekm=(k-0.5)/N + m;
      z = std::complex<double>(-PI*N*pow(ekm - q,2) , (2*PI*N*(ekm - q/2)*p - PI*m));  
      zz = exp(z);
      sum2 += zz;
    }
    res(k) = sum2;
  }
  return res;
}


/*                                           */
/* Iterate n times using the matrix operator */
/*                                           */
void iterate(LaVectorComplex & v, LaGenMatComplex & m, int n)
{
  if(n >= 1)
  {
    LaVectorComplex vTmp(v);
    for(int i=0;i<n;++i)
    {
      Blas_Mat_Vec_Mult(m,v,vTmp);
      v=vTmp;
    }
  }
}

/*                                                   */
/* Iterate n times the operator on the given vector, */
/* (represented in the eigen basis), using the       */
/* corresponding eigenvalues.                        */
/* It modifies the initial vector.                   */
/*                                                   */

void iterate(LaGenMatComplex & v, LaVectorComplex & eigVal, LaGenMatComplex & eigVecL, LaGenMatComplex & eigVecR, int n)
{
  if(n>=1)
  {
    if(v.cols() == 1 && v.rows() == eigVal.size())
      for(unsigned int i=0;i<eigVal.size();++i)
      {
        // cannot avoid the std::complex<double> for the pow function
        std::complex<double> lambda(eigVal(i).r, eigVal(i).i);
//        LaComplex scal = complexScalarProd(eigVecL.col(i),eigVecR.col(i));
        
        lambda = pow(lambda,n);
        LaComplex lambdan (lambda);
        LaComplex comp(v(i,0));
        v(i,0) = comp*lambdan*complexScalarProd(eigVecL.col(i),eigVecR.col(i));
      }
    else
      std::cerr << "Error in iterate(LaVectorComplex & v, LaVectorComplex & eigVal, int n) : v and eigVal don't have the same size." << std::endl;
  }
}

/*                                                   */
/* The same, but if the vector is a matrix (M,1)     */
/*                                                   */

void iterateWithLambda(LaGenMatComplex & v, LaVectorComplex & eigVal, LaGenMatComplex & eigVecL, LaGenMatComplex & eigVecR, double n)
{
  std::cout << "IterateWithLambda..." << std::endl;
  if(v.cols() == 1 && v.rows() == eigVal.size())
  {
    LaGenMatComplex vtmp(v);
    for(unsigned j=0;j<eigVal.size();++j)
    {
//      std::cerr << j << std::endl;
      LaComplex res(0,0);
      for(unsigned int i=0;i<eigVal.size();++i)
      {
        std::complex<double> lambda(eigVal(i).r, eigVal(i).i);
        lambda = pow(lambda,n);
        LaComplex lambdan (lambda);
        LaComplex comp(vtmp(i,0));
        LaComplex scal(eigVecR(j,i));
        res = res + comp*lambdan*scal;
//        std::cout << res << std::endl;
//        std::cerr << "---------" << std::endl;
//        std::cerr << "lambda : " << lambda << std::endl;
//        std::cerr << "lambdan : " << lambdan << std::endl;
//        std::cerr << "comp : " << comp << std::endl;
//        std::cerr <<  "scal : " << scal << std::endl;
      }
      v(j,0) = res;
    }
  }
  else
    std::cerr << "Error in iterate(LaMatComplex & v, LaVectorComplex & eigVal, int n) : v has not one column or v and eigVal don't have the same size." << std::endl;
}




void storeResult( const char * fileName, 
                  std::vector<std::vector<Particle> > mesh,
                  LaVectorDouble result,
                  unsigned int dimX,
                  unsigned int dimY,
                  std::string dataType)
{
//  LaVectorDouble result; unsigned int dimX; unsigned int dimY;
  std::fstream file(fileName,std::ios::out);
  
  if(file)
  {

    std::ostringstream oss;
    oss << "Husimi_average";
    std::string  title = oss.str();
    
    // For Hyplot-
    std::cout << "Storing data..."                                << std::endl;
    file      << "#DATA " << dataType                             << std::endl;
    file      << "#DIM " << dimX << " " << dimY                   << std::endl;
    file      << "#COLOR_BG 0"                                    << std::endl;
    file      << "#TITLE " << title                               << std::endl;
    file      << "#BEGIN"                                         << std::endl;

    if(dataType == "XYZ")
    {
      if(result.size() == mesh.size()*mesh[0].size())
      {
        for(unsigned int i = 0;i<mesh.size();++i)
        {
          fprintf(stdout,"%lf%%\r",(double)(i*100)/mesh.size());
          for(unsigned int j = 0;j<mesh[0].size();++j)
          {
            file << mesh[i][j] << " " << result(i*mesh[0].size() + j) << std::endl;
          }
        }
      }
     else
     {
       std::cerr << "Error : results and mesh have not the same size." << std::endl;
     }
     
    }
    else if(dataType == "XY")
    {
      for(unsigned int i = 0;i<mesh.size();++i)
      {
        fprintf(stdout,"%lf%%\r",(double)(i*100)/mesh.size());
        file << mesh[i][0].getQ() << " " << result(i) << std::endl;
      }
    }
    else
      std::cerr << "Wrong data type." << std::endl;
  }
  else
    std::cerr << "Error, output file already exist." << std::endl;
}



#endif
