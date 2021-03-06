#include "BOV.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* fonction that outputs a random value, with a
 * probability following a gaussian curve */
GLfloat random_gauss(GLfloat mu, GLfloat sigma);

/* fill coord with random coordinates following a uniform distribution */
void random_uniform_points(GLfloat coord[][2], GLsizei n,
                           GLfloat min[2], GLfloat max[2]);

/* creating random points following a gaussian distribution.
 * around multiple centroid (maximum 6 centroids) which
 * are uniformly */
void random_points(GLfloat coord[][2], GLsizei n);

/* create a random polygon
 * the bigger nSmooth is, the rounder it will be  */
void random_polygon(GLfloat coord[][2], GLsizei n, int nSmooth);


typedef struct point point;
struct point
{
    double angle;
    double y;
    double x;
    int index;
};

/* this function is used by 'qsort' to sort the points with respect to
 * their y-coordinate
 */
int compAngle(const void* a, const void* b);

/* this function is used by 'qsort' to sort the points with respect to
 * the angle they form with the horizontal line passing through the lowest point
 */
int compY(const void* a, const void* b);

/* returns the distance between two points
 */
double dist(int a, int b, GLfloat coord[][2]);