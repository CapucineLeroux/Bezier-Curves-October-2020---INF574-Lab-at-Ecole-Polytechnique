#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>

using namespace Eigen;
using namespace std;

/**
 * A class for dealin with the computation of Bezier curves and their properties
 **/
class Bezier
{

public:
  /**
 * An iterative implementation of the De Casteljau algorithm for computing a Bezier curve
 * 
 * @param V  the vertices of the control polygon
 * @param t   an input parameter in [0, 1]
 * 
 * @return    the point B(t), obtaining by evaluating the curve for the value 't' 
 **/
  MatrixXd de_casteljau(const MatrixXd &V, double t)
  {
    int nV = V.rows();         // number of points in the control polygon
    int degree = V.rows() - 1; // degree of the curve

    MatrixXd B = V;

    for (int r=1 ; r<=degree ; r++){
        for (int i=0 ; i<=(degree-r) ; i++){
            B.row(i) = (1.-t)*B.row(i) + t*B.row(i+1);
        }
    }

    return B.row(0);

  }

  /**
	 * Plot the curve, for t=0, dt, 2*dt, ..., 1, with a given 'resolution' <br>
   * where dt=1/(resolution-1)
   * 
   * @param resolution  number of points to be evaluated on the curve
	 */
  MatrixXd plot_curve(const MatrixXd &V, int resolution)
  {
    double dt = 1. / (resolution - 1);
    MatrixXd Bezier(resolution,3);

    for (int i=0 ; i<resolution ; i++){
        Bezier.row(i) = de_casteljau(V,i*dt);;
    }

    return Bezier;

  }

  /**
	 * Perform the subdivision (once) of the Bezier curve (with parameter t) <br>
	 * Return two Bezier curves (with 'n' control points each)
	 */
  vector<MatrixXd> subdivide(const MatrixXd &V, double t)
  {
    vector<MatrixXd> curves{}; // the result: store the 2 curves obtained after subdivision


    int n = V.rows();         // number of points in the control polygon
    int degree = V.rows() - 1; // degree of the curve

    //create 2 matrices V1 and V2 with the new smaller control polygons
    MatrixXd V1(n,3);
    MatrixXd V2(n,3);

    //Matrix to calculate the points of V1 and V2
    MatrixXd B = V;
    V1.row(0) = B.row(0);
    V2.row(degree-0) = B.row(degree);

    for (int r=1 ; r<=degree ; r++){
        for (int i=0 ; i<=(degree-r) ; i++){
            B.row(i) = (1.-t)*B.row(i) + t*B.row(i+1);
        }
        V1.row(r) = B.row(0);
        V2.row(degree-r) = B.row(degree-r);
    }

    curves.push_back(V1);
    curves.push_back(V2);

    return curves;
  }

  /**
	 * Plot the curve using recursive subdivision <br>
   * 
   * @param levels  number of levels of subdivisions
   * @return  a polyline representing the curve to be rendered: this is obtained by concatenation of
   * the control points of all subdivided curves
	 */
 MatrixXd subdivision_plot(const MatrixXd &V, int levels)
  {
    std::cout << "computing recursive subdivision " << std::endl;
    vector<MatrixXd> curves{};

    //calculates the two new control polygons, splitting V in two
    vector<MatrixXd> SV = subdivide(V,0.5);
    MatrixXd V1 = SV[0];
    MatrixXd V2 = SV[1];

    int k = 4;

    //if the max number of levels is not reached, we apply the same algorithm on both halves
    if (levels <= k){
        curves.push_back(subdivision_plot(V1,levels+1));
        curves.push_back(subdivision_plot(V2,levels+1));
    }

    //if the max number of levels is reached, we keep the two polygons as they are
    else {
        curves.push_back(V1);
        curves.push_back(V2);
    }

    //concatenates in a matrix all the control points stored in curves
    int n1 = curves[0].rows();
    int n2 = curves[1].rows();
    MatrixXd CV(n1+n2,3);
    for (int i=0 ; i<n1 ; i++){
        CV.row(i) = curves[0].row(i);
    }
    for (int i=0 ; i<n2 ; i++){

        CV.row(n1+i) = curves[1].row(i);
    }

    return CV;

  }

  /**
 * Compute the tangent of a given curve c(t) for a given parameter t0
 * 
 * @param V  the vertices of the control polygon
 * @param t0   an input parameter in [0, 1]
 * 
 * @return    the tangent at c(t0)
 **/
  MatrixXd compute_tangent(const MatrixXd &V, double t0)
  {
    int n = V.rows()-1; // number of points in the control polygon

    //initialization of the elements of the sum for i=0
    MatrixXd Db = n*pow(1.-t0,n-1)*(V.row(1)-V.row(0));
    double coeff_bin = 1.;

    //calculates the sum formula using Bernstein polynomials
    for (int i=1 ; i<n ; i++){
        coeff_bin = coeff_bin*(n-i)/i;
        Db = Db + n*coeff_bin*pow(t0,i)*pow(1.-t0,n-1-i)*(V.row(i+1)-V.row(i));
    }

    //normalize Db
    double N = Db.norm();
    Db = Db/N;

    return Db;
  }

  /**
 * Compute the normal vector of a given curve c(t) for a given parameter t0
 * 
 * @param V  the vertices of the control polygon
 * @param t0   an input parameter in [0, 1]
 * 
 * @return    the normal at c(t0)
 **/
  MatrixXd compute_normal(const MatrixXd &V, double t0)
  {
    int n = V.rows(); // number of points in the control polygon

    //calculates the derivative of the tangent in t0
    double delta_t = 0.01;
    double t1 = std::max(t0-delta_t,0.);
    double t2 = std::min(t0+delta_t,1.);

    MatrixXd Db1 = compute_tangent(V,t1);
    MatrixXd Db2 = compute_tangent(V,t2);
    MatrixXd Nb = Db2 - Db1;

    //normalizes the normal vector
    double N = Nb.norm();
    Nb = Nb/N;

    return Nb;

  }

  //computes the orthogonal vector of the normal and the tangent vector
  MatrixXd compute_v(const MatrixXd &V, double t0)
  {
      MatrixXd T0 = compute_tangent(V,t0);
      MatrixXd N0 = compute_normal(V,t0);
      Vector3d T = T0.row(0);
      Vector3d N = N0.row(0);
      Vector3d v = N.cross(T);
      MatrixXd v0(1,3);
      v0.row(0) = v;

      double n = v0.norm();
      v0 = v0/n;
      return v0;
  }

  /**
 * Compute a loop of points around a curve c(t) for a given parameter t0
 * The points belongs on a circle lying the hyperplane passing through c(t0) and orthogonal to tangent(t0)
 * 
 * @param V  the vertices of the control polygon
 * @param t0   an input parameter in [0, 1]
 * 
 * @return    a loop of vertices on the hyperplane passing through c(t0) and orthogonal to tangent(t0)
 **/
  MatrixXd compute_loop_of_vertices(const MatrixXd &V, double t0, int k, double radius)
  {

    //defines hyperplane using the normal and the other orthogonal vector of the tangent
    MatrixXd X = compute_normal(V,t0);
    MatrixXd Y = compute_v(V,t0);
    //c(t0)
    MatrixXd P = de_casteljau(V,t0);

    double theta = 2.*3.14/(k-1);
    MatrixXd loop(k,3);
    double theta_i;

    for (int i=0 ; i<k ; i++){
        theta_i = i*theta;
        loop.row(i) = radius*std::cos(theta_i)*X + radius*std::sin(theta_i)*Y + P;
    }

    return loop;
  }

};
