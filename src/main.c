#include "inputs.h"
#include <time.h>
#include <math.h>
#include <stdio.h>
#include "myRobustPredicates.h"

#ifndef TRUE
#define TRUE 1
#endif // !TRUE

double x_low, y_low;
GLfloat(*coord)[2];
bov_points_t* coordDraw;
GLsizei nPoints;
int nPoints_scan;
int nPoints_max = 500;
int nPoints_supp;
int ready;
int points_recorded;

int lowest(GLfloat coord[][2], GLsizei nPoints);

/* fills the vector 'sorted_coord' with the indices of the points ordered with respect to
 * the angle they form with the horizontal line passing through the lowest point
 */
void sort_angles(GLfloat coord[][2], GLsizei nPoints, GLuint* sorted_coord) {

	int low_point = lowest(coord, nPoints);
	sorted_coord[0] = low_point;
	x_low = coord[low_point][0];
	y_low = coord[low_point][1];

	point* arrayPoints = malloc(nPoints * sizeof(point));
	if (!arrayPoints)
		exit(0);

	for (int i = 0; i < nPoints; i++) {
		arrayPoints[i].index = i;
		arrayPoints[i].x = coord[i][0];
		arrayPoints[i].y = coord[i][1];
		if (i == low_point)
			arrayPoints[i].angle = -1000.0;
		else if (coord[i][0] < coord[low_point][0])
			arrayPoints[i].angle = M_PI / 2 - atan((coord[i][0] - coord[low_point][0]) / (coord[i][1] - coord[low_point][1]));
		else
			arrayPoints[i].angle = atan((coord[i][1] - coord[low_point][1]) / (coord[i][0] - coord[low_point][0]));
	}

	qsort(&arrayPoints[0], nPoints, sizeof(arrayPoints[0]), compAngle);

	int k = 0;
	for (int i = 0; i < nPoints; i++) {
		if (i != nPoints - 1 && arrayPoints[i].angle == arrayPoints[i + 1].angle)
			continue;
		sorted_coord[k] = arrayPoints[i].index;
		k++;
	}

	free(arrayPoints);

	return;
}

/* fills the vector 'sorted_coord' with the indices of the points ordered with respect to
 * their y-coordinate
 */
void sort_y(GLfloat coord[][2], GLsizei nPoints, GLuint* sorted_coord) {
	point* arrayPoints = malloc(nPoints * sizeof(point));
	if (!arrayPoints)
		exit(0);

	for (int i = 0; i < nPoints; i++) {
		arrayPoints[i].index = i;
		arrayPoints[i].y = coord[i][1];
	}

	qsort(&arrayPoints[0], nPoints, sizeof(arrayPoints[0]), compY);

	for (int i = 0; i < nPoints; i++)
		sorted_coord[i] = arrayPoints[i].index;

	free(arrayPoints);

	return;
}

/* this function is used by 'qsort' to sort the points with respect to
 *the angle they form with the horizontal line passing through the lowest point
 */
int compAngle(const void* a, const void* b) {
	double diff_angle, diff_dist;
	diff_angle = ((point*)a)->angle - ((point*)b)->angle;
	if (diff_angle > 0.0)
		return 1;
	else if (diff_angle < 0.0)
		return -1;
	else {
		diff_dist = sqrt(pow(((point*)a)->x - x_low, 2) + pow(((point*)a)->y - y_low, 2))
			- sqrt(pow(((point*)b)->x - x_low, 2) + pow(((point*)b)->y - y_low, 2));
		if (diff_dist < 0.0)
			return 1;
		else
			return -1;
	}
}

/* this function is used by 'qsort' to sort the points with respect to
 * their y-coordinate
 */
int compY(const void* a, const void* b) {
	double diff = ((point*)a)->y - ((point*)b)->y;
	return diff >= 0.0 ? 1 : -1;
}

/* takes the coordinates of all points and 3 indices as inputs
 * returns 1 if the points 'a','b','c' define a counterclockwise oriented surface
 *        -1 if the points 'a','b','c' define a clockwise oriented surface
 *         0 if the points are colinear
 */
int orientation(int a, int b, int c, GLfloat coord[][2]) {
	double coord_a[2] = { (double)coord[a][0], (double)coord[a][1] };
	double coord_b[2] = { (double)coord[b][0], (double)coord[b][1] };
	double coord_c[2] = { (double)coord[c][0], (double)coord[c][1] };
	double s = my_orient2d(coord_a, coord_b, coord_c);
	return s > 0.0 ? 1 : s < 0.0 ? -1 : 0;
}

/* searchs in the array 'ch' until the 'count'-th element
 * returns 1 if 'number' is in the array
 *         0 if 'number' is not in the array
 */
int isInCH(int number, int* ch, int count) {
	for (int p = 0; p < count; p++) {
		if (ch[p] == number)
			return 1;
	}
	return 0;
}

/* searchs the point with the smallest y-component
 * among a set of 'nPoints' points
 */
int lowest(GLfloat coord[][2], GLsizei nPoints) {
	int low_point = 0;
	for (int i = 0; i < nPoints; i++) {
		if (coord[i][1] < coord[low_point][1])
			low_point = i;
	}
	return low_point;
}

/* returns the number of points which belong to the Convex Hull
 */
int goodSize(int* ch, GLsizei nPoints) {
	int n_ch = 0;
	for (int i = 1; i < nPoints; i++) {
		if (ch[i] < 0) {
			n_ch = i;
			break;
		}
	}
	n_ch++;

	return n_ch;
}

/* fills the vector 'ch_goodSize' with the 'n_ch' first elements of 'ch'
 */
void fill(int* ch, int n_ch, int* ch_goodSize) {
	for (int i = 0; i < n_ch - 1; i++)
		ch_goodSize[i] = ch[i];
	ch_goodSize[n_ch - 1] = ch_goodSize[0];

	return;
}

/* searchs the Convex Hull corresponding to a set of 'nPoints' points
 * thanks to the "Jarvis March" algorithm
 */
void convex_hull_jarvis_march(GLfloat coord[][2], GLsizei nPoints, int* ch) {

	printf("\n---------- JARVIS MARCH ALGORITHM ----------\n");

	int ref = 0, next = 0, count = 0;

	// The first 'ref' point is the lowest one
	ref = lowest(coord, nPoints);

	while (TRUE) {
		// The 'ref' point is added to the Convex Hull
		ch[count] = ref;
		count++;

		// The 'next' candidate is chosen: by default, we take the point following the 'ref' one
		// and we exclude all the points that belong to the Convex Hull
		if (ref == nPoints - 1)
			next = 0;
		else
			next = ref + 1;
		while (isInCH(next, ch, count) && (next != nPoints - 1))
			next++;

		// Loop on all the points to search a better candidate than 'next'
		for (int i = 0; i < nPoints; i++) {
			// We exclude:
			//    - the points that already belong to the Convex Hull except the initial one (to close the loop)
			//    - the reference point
			//    - the current candidate 'next'
			if ((isInCH(i, ch, count) && ch[0] != i) || i == next || i == ref)
				continue;

			// If the surface formed by the 3 points (ref, i and next) is
			//    - counterclockwise: 'i' is a better candidate than 'next' so we replace this latter
			//    - clockwise: 'next' is still the better candidate, we continue the iteration
			if (orientation(ref, i, next, coord) > 0) {
				next = i;
				continue;
			}
		}

		// At the end of the loop, we found the 'next' point which will belong to
		// the Convex Hull. This will be the new 'ref' point
		ref = next;

		// If the 'ref' point is the lowest one, we closed the Convex Hull
		if (ref == ch[0])
			break;
	}

	return;
}

/* searchs the Convex Hull corresponding to a set of 'nPoints' points
 * thanks to the "Graham Scan" algorithm by sorting the points with respect to
 * the angle they form with the horizontal line passing through the lowest point
 */
void convex_hull_graham_scan_angle(GLfloat coord[][2], GLsizei nPoints, int* ch) {

	printf("\n---------- GRAHAM SCAN ALGORITHM (sort by angle) ----------\n");

	GLuint* sorted_coord = malloc(nPoints * sizeof(GLuint));
	if (!sorted_coord)
		exit(0);

	// First, we sort the point with respect to their angle
	clock_t t1_sort_angle = clock();
	sort_angles(coord, nPoints, sorted_coord);
	clock_t t2_sort_angle = clock();
	printf("TOTAL SORTING TIME (angle): %f [s]\n", (float)(t2_sort_angle - t1_sort_angle) / CLOCKS_PER_SEC);

	int* stack = malloc(nPoints * sizeof(int));
	if (!stack)
		exit(0);
	int top;
	for (top = 0; top < 3; top++)
		stack[top] = sorted_coord[top];
	top--;

	for (int i = 3; i < nPoints; i++) {
		if (sorted_coord[i] < 0)
			break;

		while (orientation(stack[top - 1], stack[top], sorted_coord[i], coord) <= 0) {
			stack[top] = -1;
			top--;
		}
		top++;
		stack[top] = sorted_coord[i];
	}

	// We put the stack in the half Convex Hull
	for (int j = 0; j < nPoints; j++)
		ch[j] = stack[j];

	free(sorted_coord);
	free(stack);

	return;
}

/* searchs the Convex Hull corresponding to a set of 'nPoints' points
 * thanks to the "Graham Scan" algorithm by sorting the points with respect to
 * their y-coordinate
 */
void convex_hull_graham_scan_y(GLfloat coord[][2], GLsizei nPoints, int* ch) {

	printf("\n---------- GRAHAM SCAN ALGORITHM (sort by y) ----------\n");

	int top, i;
	int* ch_r = malloc(nPoints * sizeof(int));
	if (!ch_r)
		exit(0);
	int* ch_l = malloc(nPoints * sizeof(int));
	if (!ch_l)
		exit(0);
	GLuint* sorted_coord = malloc(nPoints * sizeof(GLuint));
	if (!sorted_coord)
		exit(0);

	// First, we sort the point with respect to their y-coordinate
	clock_t t1_sort_y = clock();
	sort_y(coord, nPoints, sorted_coord);
	clock_t t2_sort_y = clock();
	printf("TOTAL SORTING TIME (y): %f [s]\n", (float)(t2_sort_y - t1_sort_y) / CLOCKS_PER_SEC);

	/* -------------------- RIGHT PART -------------------- */
	// We find the half Convex Hull from the lowest to the highest point
	int* stack_r = malloc(nPoints * sizeof(int));
	if (!stack_r)
		exit(0);

	// We put the 3 bottom points in the stack
	for (top = 0; top < 3; top++)
		stack_r[top] = sorted_coord[top];
	top--;

	// While the 3 points of the stack are clockwise-oriented, we remove the points
	i = 3;
	while (orientation(stack_r[0], stack_r[1], stack_r[2], coord) < 0) {
		stack_r[1] = stack_r[2];
		stack_r[2] = sorted_coord[i];
		i++;
	}

	// From bottom to top: if the next candidate forms a counterclockwise turn, we keep it
	// else, we remove the wrong point from the stack
	for (; i < nPoints; i++) {
		while (orientation(stack_r[top - 1], stack_r[top], sorted_coord[i], coord) <= 0) {
			if (top == 1) {
				stack_r[top] = sorted_coord[i];
				i++;
			}
			else {
				stack_r[top] = -1;
				top--;
			}
		}
		top++;
		stack_r[top] = sorted_coord[i];
	}

	// We put the stack in the half Convex Hull
	for (int j = 0; j < nPoints; j++)
		ch_r[j] = stack_r[j];
	/* -------------------- RIGHT PART -------------------- */

	/* -------------------- LEFT  PART -------------------- */
	// We find the half Convex Hull from the highest to the lowest point
	int* stack_l = malloc(nPoints * sizeof(int));
	if (!stack_l)
		exit(0);

	// We put the 3 top points in the stack
	for (top = 0; top < 3; top++)
		stack_l[top] = sorted_coord[nPoints - 1 - top];
	top--;

	// While the 3 points of the stack are clockwise-oriented, we remove the points
	i = nPoints - 4;
	while (orientation(stack_l[0], stack_l[1], stack_l[2], coord) < 0) {
		stack_l[1] = stack_l[2];
		stack_l[2] = sorted_coord[i];
		i--;
	}

	// From top to bottom: if the next candidate forms a counterclockwise turn, we keep it
	// else, we remove the wrong point from the stack
	for (; i >= 0; i--) {
		while (orientation(stack_l[top - 1], stack_l[top], sorted_coord[i], coord) <= 0) {
			if (top == 1) {
				stack_l[top] = sorted_coord[i];
				i--;
			}
			else {
				stack_l[top] = -1;
				top--;
			}
		}
		top++;
		stack_l[top] = sorted_coord[i];
	}

	// We put the stack in the half Convex Hull
	for (int j = 0; j < nPoints; j++)
		ch_l[j] = stack_l[j];
	/* -------------------- LEFT  PART -------------------- */

	/* ---------------------- MERGING --------------------- */
	int fr, fl;
	// First, we put the half Convex Hull from right into the vector 'ch'
	for (fr = 0; fr < nPoints; fr++) {
		if (ch_r[fr] < 0)
			break;
		else
			ch[fr] = ch_r[fr];
	}

	// Then, the half Convex Hull from left is added following the first half
	fr--;
	for (fl = 1; fr + fl - 1 < nPoints; fl++) {
		if (ch_l[fl] < 0)
			break;
		else
			ch[fr + fl - 1] = ch_l[fl - 1];
	}
	/* ---------------------- MERGING --------------------- */

	free(sorted_coord);
	free(ch_r); free(ch_l);
	free(stack_r); free(stack_l);

	return;
}

/* same algorithm than 'convex_hull_jarvis_march' with animations
 */
void convex_hull_jarvis_march_animation(GLfloat coord[][2], GLsizei nPoints, int* ch_index, bov_window_t* window) {

	printf("\n---------- JARVIS MARCH ALGORITHM ----------\n");
	int ref = 0, next = 0, count = 0;
	clock_t t;

	// definition of the 'pause time' as a function of nPoints
	int pause;
	if (nPoints <= 25)
		pause = 200;
	else if (nPoints > 25 && nPoints <= 50)
		pause = 100;
	else if (nPoints > 50 && nPoints <= 200)
		pause = 40;
	else if (nPoints > 200)
		pause = 10;

	// we search the lowest point to display the text below
	ref = lowest(coord, nPoints);
	bov_text_t* text_algo = bov_text_new((GLubyte[]) { "First, we find the lowest point:\n\tit belongs to the Convex Hull" }, GL_DYNAMIC_DRAW);
	bov_text_set_fontsize(text_algo, 0.1);
	bov_text_set_boldness(text_algo, 0.35);
	bov_text_set_color(text_algo, (GLfloat[4]) { 0.9, 0.2, 0.1, 1.0 });
	bov_text_set_outline_width(text_algo, 0.8);
	bov_text_set_outline_color(text_algo, (GLfloat[4]) { 0.0, 0.0, 0.0, 1.0 });
	bov_text_set_pos(text_algo, (GLfloat[2]) { -0.8, coord[ref][1] - 0.15 });
	ref = 0;

	// we find the lowest point with animations
	bov_points_t* current_ref_point = bov_points_new(coord[0], 1, GL_DYNAMIC_DRAW);
	bov_points_set_color(current_ref_point, (GLfloat[4]) { 1.0, 0.0, 0.0, 1.0 });
	bov_points_set_width(current_ref_point, 0.020);
	bov_points_t* current_node_point = bov_points_new(coord[0], 1, GL_DYNAMIC_DRAW);
	bov_points_set_color(current_node_point, (GLfloat[4]) { 0.0, 1.0, 0.0, 1.0 });
	bov_points_set_width(current_node_point, 0.020);
	
	int i = 0;
	while (!bov_window_should_close(window) && i < nPoints) {
		bov_points_update(current_node_point, coord[i], 1);
		t = clock();
		
		while (!bov_window_should_close(window) && difftime(clock(), t) < pause) {
			bov_text_draw(window, text_algo);
			bov_points_draw(window, coordDraw, 0, nPoints);
			bov_points_draw(window, current_node_point, 0, 1);
			bov_points_draw(window, current_ref_point, 0, 1);
			bov_window_update(window);
		}
		
		if (coord[i][1] < coord[ref][1]) {
			ref = i;
			bov_points_update(current_ref_point, coord[ref], 1);
		}
		i++;
	}
	bov_points_update(current_ref_point, coord[ref], 1);

	GLfloat(*ch_coord)[2] = malloc(nPoints * sizeof(coord[0]));
	if (!ch_coord)
		exit(0);
	GLfloat(*temp)[2] = malloc(2 * sizeof(coord[0]));
	if (!temp)
		exit(0);
	temp[0][0] = coord[ref][0];
	temp[0][1] = coord[ref][1];
	if (ref == 0) {
		temp[1][0] = coord[1][0];
		temp[1][1] = coord[1][1];
	}
	else {
		temp[1][0] = coord[0][0];
		temp[1][1] = coord[0][1];
	}
	
	// we update the text and we can make a little pause
	bov_text_update(text_algo, (GLubyte[]) { "Now, we can iterate on each point" });
	t = clock();
	while (!bov_window_should_close(window) && difftime(clock(), t) < 400) {
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, current_ref_point, 0, nPoints);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}

	bov_points_t* chDraw = bov_points_new(ch_coord, count + 1, GL_DYNAMIC_DRAW);
	bov_points_set_color(chDraw, (GLfloat[4]) { 0.0, 1.0, 0.0, 1.0 });
	bov_points_set_width(chDraw, 0.010);
	bov_points_t* tempDraw = bov_points_new(temp, 2, GL_DYNAMIC_DRAW);
	bov_points_set_color(tempDraw, (GLfloat[4]) { 1.0, 0.078, 0.577, 1.0 });
	bov_points_set_width(tempDraw, 0.010);
	
	// we iterate on the points until the convex hull is found
	while (TRUE) {
		ch_index[count] = ref;
		ch_coord[count][0] = coord[ref][0];
		ch_coord[count][1] = coord[ref][1];
		count++;
		bov_points_update(chDraw, ch_coord, count);

		if (ref == nPoints - 1)
			next = 0;
		else
			next = ref + 1;
		while (isInCH(next, ch_index, count) && (next != nPoints - 1))
			next++;

		temp[0][0] = coord[ref][0];
		temp[0][1] = coord[ref][1];
		temp[1][0] = coord[next][0];
		temp[1][1] = coord[next][1];
		bov_points_update(tempDraw, temp, 2);

		for (int i = 0; i < nPoints; i++) {

			if ((isInCH(i, ch_index, count) && ch_index[0] != i) || i == next || i == ref)
				continue;

			ch_coord[count][0] = coord[next][0];
			ch_coord[count][1] = coord[next][1];
			bov_points_update(chDraw, ch_coord, count + 1);

			temp[0][0] = coord[i][0];
			temp[0][1] = coord[i][1];
			temp[1][0] = coord[next][0];
			temp[1][1] = coord[next][1];
			bov_points_update(tempDraw, temp, 2);
			t = clock();
			while (!bov_window_should_close(window) && difftime(clock(), t) < pause) {
				bov_points_draw(window, coordDraw, 0, nPoints);
				bov_points_draw(window, current_ref_point, 0, nPoints);
				bov_line_strip_draw(window, tempDraw, 0, 2);
				bov_line_strip_draw(window, chDraw, 0, count + 1);
				bov_text_draw(window, text_algo);
				bov_window_update(window);
			}

			if (orientation(ref, i, next, coord) > 0) {
				next = i;
				continue;
			}
		}

		ref = next;
		if (ref == ch_index[0])
			break;
	}

	bov_text_update(text_algo, (GLubyte[]) { "Convex Hull found !" });
	printf("----- DONE -----");
	while (!bov_window_should_close(window)) {
		bov_line_loop_draw(window, chDraw, 0, nPoints);
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}

	bov_text_delete(text_algo);
	bov_points_delete(current_ref_point);
	bov_points_delete(current_node_point);
	bov_points_delete(chDraw);
	bov_points_delete(tempDraw);
	free(ch_coord);
	free(temp);

	return;
}

/* same algorithm than 'convex_hull_graham_scan_angle' with animations
 */
void convex_hull_graham_scan_angle_animation(GLfloat coord[][2], GLsizei nPoints, int* ch_index, bov_window_t* window) {

	printf("\n---------- GRAHAM SCAN ALGORITHM (sort by angle) ----------\n");

	GLuint* sorted_coord = malloc(nPoints * sizeof(GLuint));
	if (!sorted_coord)
		exit(0);
	int ref = 0, count = 0;
	clock_t t;

	// definition of the 'pause time' as a function of nPoints
	int pause;
	if (nPoints <= 25)
		pause = 200;
	else if (nPoints > 25 && nPoints <= 50)
		pause = 100;
	else if (nPoints > 50 && nPoints <= 200)
		pause = 40;
	else if (nPoints > 200)
		pause = 10;

	// we search the lowest point to display the text below
	ref = lowest(coord, nPoints);
	bov_text_t* text_algo = bov_text_new((GLubyte[]) { "First, we find the lowest point:\n\tit belongs to the Convex Hull" }, GL_DYNAMIC_DRAW);
	bov_text_set_fontsize(text_algo, 0.1);
	bov_text_set_boldness(text_algo, 0.35);
	bov_text_set_color(text_algo, (GLfloat[4]) { 0.9, 0.2, 0.1, 1.0 });
	bov_text_set_outline_width(text_algo, 0.8);
	bov_text_set_outline_color(text_algo, (GLfloat[4]) { 0.0, 0.0, 0.0, 1.0 });
	bov_text_set_pos(text_algo, (GLfloat[2]) { -0.8, coord[ref][1] - 0.15 });
	ref = 0;

	// we find the lowest point with animations
	bov_points_t* current_ref_point = bov_points_new(coord[0], 1, GL_DYNAMIC_DRAW);
	bov_points_set_color(current_ref_point, (GLfloat[4]) { 1.0, 0.0, 0.0, 1.0 });
	bov_points_set_width(current_ref_point, 0.020);
	bov_points_t* current_node_point = bov_points_new(coord[0], 1, GL_DYNAMIC_DRAW);
	bov_points_set_color(current_node_point, (GLfloat[4]) { 0.0, 1.0, 0.0, 1.0 });
	bov_points_set_width(current_node_point, 0.020);

	int i = 0;
	while (!bov_window_should_close(window) && i < nPoints) {
		bov_points_update(current_node_point, coord[i], 1);
		t = clock();
		while (!bov_window_should_close(window) && difftime(clock(), t) < pause) {
			bov_text_draw(window, text_algo);
			bov_points_draw(window, coordDraw, 0, nPoints);
			bov_points_draw(window, current_node_point, 0, 1);
			bov_points_draw(window, current_ref_point, 0, 1);
			bov_window_update(window);
		}

		if (coord[i][1] < coord[ref][1]) {
			ref = i;
			bov_points_update(current_ref_point, coord[ref], 1);
		}
		i++;
	}
	bov_points_update(current_ref_point, coord[ref], 1);

	bov_text_update(text_algo, (GLubyte[]) { "Here is the actual ordering of the points" });
	t = clock();
	while (!bov_window_should_close(window) && difftime(clock(), t) < 500) {
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, current_ref_point, 0, 1);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}

	for (int i = 0; i < nPoints; i++) {
		t = clock();
		while (!bov_window_should_close(window) && difftime(clock(), t) < pause) {
			bov_line_strip_draw(window, coordDraw, 0, i);
			bov_points_draw(window, coordDraw, 0, nPoints);
			bov_points_draw(window, current_ref_point, 0, 1);
			bov_text_draw(window, text_algo);
			bov_window_update(window);
		}
	}

	bov_text_update(text_algo, (GLubyte[]) { "This is not optimal ...\n\tLet's sort them following 'theta' : 3" });
	t = clock();
	while (!bov_window_should_close(window) && difftime(clock(), t) < 1000) {
		bov_line_strip_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, current_ref_point, 0, 1);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}

	bov_text_update(text_algo, (GLubyte[]) { "This is not optimal ...\n\tLet's sort them following 'theta' : 2" });
	t = clock();
	while (!bov_window_should_close(window) && difftime(clock(), t) < 1000) {
		bov_line_strip_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, current_ref_point, 0, 1);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}

	bov_text_update(text_algo, (GLubyte[]) { "This is not optimal ...\n\tLet's sort them following 'theta' : 1" });
	t = clock();
	while (!bov_window_should_close(window) && difftime(clock(), t) < 1000) {
		bov_line_strip_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, current_ref_point, 0, 1);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}

	sort_angles(coord, nPoints, sorted_coord);
	bov_order_t* order_sort = bov_order_new(sorted_coord, nPoints, GL_DYNAMIC_DRAW);
	GLfloat(*closeLoopCoord)[2] = malloc(2 * sizeof(coord[0]));
	if (!closeLoopCoord)
		exit(0);
	closeLoopCoord[0][0] = coord[sorted_coord[nPoints - 1]][0];
	closeLoopCoord[0][1] = coord[sorted_coord[nPoints - 1]][1];
	closeLoopCoord[1][0] = coord[sorted_coord[0]][0];
	closeLoopCoord[1][1] = coord[sorted_coord[0]][1];
	bov_points_t* closeLoopDraw = bov_points_new(closeLoopCoord, 2, GL_DYNAMIC_DRAW);
	bov_points_set_color(closeLoopDraw, (GLfloat[4]) { 0.0, 0.0, 0.0, 1.0 });
	bov_points_set_width(closeLoopDraw, 0.005);
	bov_text_update(text_algo, (GLubyte[]) { "This is much better ! We can do the iteration" });
	t = clock();
	while (!bov_window_should_close(window) && difftime(clock(), t) < 2000) {
		bov_line_strip_draw_with_order(window, coordDraw, order_sort, 0, nPoints);
		bov_line_strip_draw(window, closeLoopDraw, 0, 2);
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, current_ref_point, 0, 1);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}

	GLfloat(*ch_coord)[2] = malloc(nPoints * sizeof(coord[0]));
	if (!ch_coord)
		exit(0);
	GLfloat(*temp)[2] = malloc(2 * sizeof(coord[0]));
	if (!temp)
		exit(0);
	bov_points_t* chDraw = bov_points_new(ch_coord, count + 1, GL_DYNAMIC_DRAW);
	bov_points_set_color(chDraw, (GLfloat[4]) { 0.0, 1.0, 0.0, 1.0 });
	bov_points_set_width(chDraw, 0.010);
	bov_points_t* tempDraw = bov_points_new(temp, 2, GL_DYNAMIC_DRAW);
	bov_points_set_color(tempDraw, (GLfloat[4]) { 1.0, 0.078, 0.577, 1.0 });
	bov_points_set_width(tempDraw, 0.010);

	int* stack = malloc(nPoints * sizeof(int));
	if (!stack)
		exit(0);
	int top;
	for (top = 0; top < 3; top++) {
		stack[top] = sorted_coord[top];
		ch_coord[top][0] = coord[sorted_coord[top]][0];
		ch_coord[top][1] = coord[sorted_coord[top]][1];
	}
	top--;
	bov_points_update(chDraw, ch_coord, top);

	temp[0][0] = ch_coord[1][0];
	temp[0][1] = ch_coord[1][1];
	temp[1][0] = ch_coord[2][0];
	temp[1][1] = ch_coord[2][1];
	bov_points_update(tempDraw, temp, 2);

	t = clock();
	while (!bov_window_should_close(window) && difftime(clock(), t) < pause) {
		bov_line_strip_draw_with_order(window, coordDraw, order_sort, 0, nPoints);
		bov_line_strip_draw(window, closeLoopDraw, 0, 2);
		bov_line_strip_draw(window, chDraw, 0, top + 1);
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, current_ref_point, 0, 1);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}

	for (int i = 3; i < nPoints; i++) {
		while (orientation(stack[top - 1], stack[top], sorted_coord[i], coord) <= 0) {
			stack[top] = -1;
			top--;

			temp[0][0] = ch_coord[top][0];
			temp[0][1] = ch_coord[top][1];
			temp[1][0] = ch_coord[top + 1][0];
			temp[1][1] = ch_coord[top + 1][1];
			bov_points_update(tempDraw, temp, 2);

			t = clock();
			while (!bov_window_should_close(window) && difftime(clock(), t) < pause) {
				bov_line_strip_draw_with_order(window, coordDraw, order_sort, 0, nPoints);
				bov_line_strip_draw(window, closeLoopDraw, 0, 2);
				bov_line_strip_draw(window, chDraw, 0, top + 1);
				bov_line_strip_draw(window, tempDraw, 0, 2);
				bov_points_draw(window, coordDraw, 0, nPoints);
				bov_points_draw(window, current_ref_point, 0, 1);
				bov_text_draw(window, text_algo);
				bov_window_update(window);
			}
		}

		top++;
		stack[top] = sorted_coord[i];
		ch_coord[top][0] = coord[sorted_coord[i]][0];
		ch_coord[top][1] = coord[sorted_coord[i]][1];
		bov_points_update(chDraw, ch_coord, top);

		t = clock();
		while (!bov_window_should_close(window) && difftime(clock(), t) < pause) {
			bov_line_strip_draw_with_order(window, coordDraw, order_sort, 0, nPoints);
			bov_line_strip_draw(window, closeLoopDraw, 0, 2);
			bov_line_strip_draw(window, chDraw, 0, top);
			bov_line_strip_draw(window, tempDraw, 0, 2);
			bov_points_draw(window, coordDraw, 0, nPoints);
			bov_points_draw(window, current_ref_point, 0, 1);
			bov_text_draw(window, text_algo);
			bov_window_update(window);
		}
	}
	bov_points_update(chDraw, ch_coord, top + 1);
	bov_text_update(text_algo, (GLubyte[]) { "Convex Hull found !" });
	printf("----- DONE -----");
	while (!bov_window_should_close(window)) {
		bov_line_loop_draw(window, chDraw, 0, nPoints);
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}

	free(sorted_coord);
	bov_text_delete(text_algo);
	bov_points_delete(current_ref_point);
	bov_points_delete(current_node_point);
	bov_points_delete(closeLoopDraw);
	bov_points_delete(chDraw);
	bov_points_delete(tempDraw);
	bov_order_delete(order_sort);
	free(closeLoopCoord);
	free(ch_coord);
	free(temp);

	return;
}

/* same algorithm than 'convex_hull_graham_scan_y' with animations
 */
void convex_hull_graham_scan_y_animation(GLfloat coord[][2], GLsizei nPoints, int* ch_index, bov_window_t* window) {

	printf("\n---------- GRAHAM SCAN ALGORITHM (sort by y) ----------\n");

	GLuint* sorted_coord = malloc(nPoints * sizeof(GLuint));
	if (!sorted_coord)
		exit(0);
	int ref = 0, count = 0, i = 0;
	clock_t t;

	// definition of the 'pause time' as a function of nPoints
	int pause;
	if (nPoints <= 25)
		pause = 200;
	else if (nPoints > 25 && nPoints <= 50)
		pause = 100;
	else if (nPoints > 50 && nPoints <= 200)
		pause = 40;
	else if (nPoints > 200)
		pause = 10;

	// we search the lowest point to display the text below
	ref = lowest(coord, nPoints);
	bov_text_t* text_algo = bov_text_new((GLubyte[]) { "First, we find the extreme points:\nthey belong to the Convex Hull" }, GL_DYNAMIC_DRAW);
	bov_text_set_fontsize(text_algo, 0.1);
	bov_text_set_boldness(text_algo, 0.35);
	bov_text_set_color(text_algo, (GLfloat[4]) { 0.9, 0.2, 0.1, 1.0 });
	bov_text_set_outline_width(text_algo, 0.8);
	bov_text_set_outline_color(text_algo, (GLfloat[4]) { 0.0, 0.0, 0.0, 1.0 });
	bov_text_set_pos(text_algo, (GLfloat[2]) { -0.8, coord[ref][1] - 0.15 });

	// we find the lowest point with animations
	bov_points_t* current_low_point = bov_points_new(coord[0], 1, GL_DYNAMIC_DRAW);
	bov_points_set_color(current_low_point, (GLfloat[4]) { 1.0, 0.0, 0.0, 1.0 });
	bov_points_set_width(current_low_point, 0.020);
	bov_points_t* current_high_point = bov_points_new(coord[0], 1, GL_DYNAMIC_DRAW);
	bov_points_set_color(current_high_point, (GLfloat[4]) { 0.0, 0.0, 1.0, 1.0 });
	bov_points_set_width(current_high_point, 0.020);
	bov_points_t* current_node_point = bov_points_new(coord[0], 1, GL_DYNAMIC_DRAW);
	bov_points_set_color(current_node_point, (GLfloat[4]) { 0.0, 1.0, 0.0, 1.0 });
	bov_points_set_width(current_node_point, 0.020);

	i = 0;
	int low = 0, high = 0;
	while (!bov_window_should_close(window) && i < nPoints) {
		bov_points_update(current_node_point, coord[i], 1);
		t = clock();

		while (!bov_window_should_close(window) && difftime(clock(), t) < pause) {
			bov_text_draw(window, text_algo);
			bov_points_draw(window, coordDraw, 0, nPoints);
			bov_points_draw(window, current_node_point, 0, 1);
			bov_points_draw(window, current_low_point, 0, 1);
			bov_points_draw(window, current_high_point, 0, 1);
			bov_window_update(window);
		}

		if (coord[i][1] < coord[low][1]) {
			low = i;
			bov_points_update(current_low_point, coord[low], 1);
		}

		if (coord[i][1] > coord[high][1]) {
			high = i;
			bov_points_update(current_high_point, coord[high], 1);
		}
		i++;
	}
	bov_points_update(current_low_point, coord[low], 1);
	bov_points_update(current_high_point, coord[high], 1);

	bov_text_update(text_algo, (GLubyte[]) { "Here is the actual ordering of the points" });
	t = clock();
	while (!bov_window_should_close(window) && difftime(clock(), t) < 500) {
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, current_low_point, 0, 1);
		bov_points_draw(window, current_high_point, 0, 1);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}

	for (i = 0; i < nPoints; i++) {
		t = clock();
		while (!bov_window_should_close(window) && difftime(clock(), t) < pause) {
			bov_line_strip_draw(window, coordDraw, 0, i);
			bov_points_draw(window, coordDraw, 0, nPoints);
			bov_points_draw(window, current_low_point, 0, 1);
			bov_points_draw(window, current_high_point, 0, 1);
			bov_text_draw(window, text_algo);
			bov_window_update(window);
		}
	}

	bov_text_update(text_algo, (GLubyte[]) { "This is not optimal ...\n\tLet's sort them following 'y' : 3" });
	t = clock();
	while (!bov_window_should_close(window) && difftime(clock(), t) < 1000) {
		bov_line_strip_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, current_low_point, 0, 1);
		bov_points_draw(window, current_high_point, 0, 1);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}

	bov_text_update(text_algo, (GLubyte[]) { "This is not optimal ...\n\tLet's sort them following 'y' : 2" });
	t = clock();
	while (!bov_window_should_close(window) && difftime(clock(), t) < 1000) {
		bov_line_strip_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, current_low_point, 0, 1);
		bov_points_draw(window, current_high_point, 0, 1);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}

	bov_text_update(text_algo, (GLubyte[]) { "This is not optimal ...\n\tLet's sort them following 'y' : 1" });
	t = clock();
	while (!bov_window_should_close(window) && difftime(clock(), t) < 1000) {
		bov_line_strip_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, current_low_point, 0, 1);
		bov_points_draw(window, current_high_point, 0, 1);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}

	sort_y(coord, nPoints, sorted_coord);
	bov_order_t* order_sort = bov_order_new(sorted_coord, nPoints, GL_DYNAMIC_DRAW);	

	bov_text_update(text_algo, (GLubyte[]) { "This is much better ! We can do the iteration" });
	t = clock();
	while (!bov_window_should_close(window) && difftime(clock(), t) < 1000) {
		bov_line_strip_draw_with_order(window, coordDraw, order_sort, 0, nPoints);
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, current_low_point, 0, 1);
		bov_points_draw(window, current_high_point, 0, 1);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}

	GLfloat(*ch_coord_l)[2] = malloc(nPoints * sizeof(coord[0]));
	if (!ch_coord_l)
		exit(0);
	GLfloat(*ch_coord_r)[2] = malloc(nPoints * sizeof(coord[0]));
	if (!ch_coord_r)
		exit(0);
	GLfloat(*temp)[2] = malloc(2 * sizeof(coord[0]));
	if (!temp)
		exit(0);
	bov_points_t* chDraw_l = bov_points_new(ch_coord_r, count + 1, GL_DYNAMIC_DRAW);
	bov_points_set_color(chDraw_l, (GLfloat[4]) { 0.0, 0.0, 1.0, 1.0 });
	bov_points_set_width(chDraw_l, 0.010);
	bov_points_t* chDraw_r = bov_points_new(ch_coord_r, count + 1, GL_DYNAMIC_DRAW);
	bov_points_set_color(chDraw_r, (GLfloat[4]) { 1.0, 0.0, 0.0, 1.0 });
	bov_points_set_width(chDraw_r, 0.010);
	bov_points_t* tempDraw = bov_points_new(temp, 2, GL_DYNAMIC_DRAW);
	bov_points_set_color(tempDraw, (GLfloat[4]) { 1.0, 0.078, 0.577, 1.0 });
	bov_points_set_width(tempDraw, 0.010);

	int top;

	/****************** RIGHT PART ******************/
	int* stack_r = malloc(nPoints * sizeof(int));
	if (!stack_r)
		exit(0);
	bov_text_update(text_algo, (GLubyte[]) { "Let's begin with the right part" });
	t = clock();
	while (!bov_window_should_close(window) && difftime(clock(), t) < 1000) {
		bov_line_strip_draw_with_order(window, coordDraw, order_sort, 0, nPoints);
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, current_low_point, 0, 1);
		bov_points_draw(window, current_high_point, 0, 1);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}

	for (top = 0; top < 3; top++) {
		stack_r[top] = sorted_coord[top];
		ch_coord_r[top][0] = coord[sorted_coord[top]][0];
		ch_coord_r[top][1] = coord[sorted_coord[top]][1];
	}
	top--;
	i = 3;
	while (orientation(stack_r[0], stack_r[1], stack_r[2], coord) < 0) {
		stack_r[1] = stack_r[2];
		ch_coord_r[1][0] = ch_coord_r[2][0];
		ch_coord_r[1][1] = ch_coord_r[2][1];
		stack_r[2] = sorted_coord[i];
		ch_coord_r[2][0] = coord[sorted_coord[i]][0];
		ch_coord_r[2][1] = coord[sorted_coord[i]][1];
		bov_points_update(chDraw_r, ch_coord_r, top);
		temp[0][0] = ch_coord_r[1][0];
		temp[0][1] = ch_coord_r[1][1];
		temp[1][0] = ch_coord_r[2][0];
		temp[1][1] = ch_coord_r[2][1];
		bov_points_update(tempDraw, temp, 2);
		t = clock();
		while (!bov_window_should_close(window) && difftime(clock(), t) < pause) {
			bov_line_strip_draw_with_order(window, coordDraw, order_sort, 0, nPoints);
			bov_line_strip_draw(window, chDraw_r, 0, nPoints);
			bov_points_draw(window, coordDraw, 0, nPoints);
			bov_points_draw(window, current_low_point, 0, 1);
			bov_points_draw(window, current_high_point, 0, 1);
			bov_text_draw(window, text_algo);
			bov_window_update(window);
		}
		i++;
	}
	
	temp[0][0] = ch_coord_r[1][0];
	temp[0][1] = ch_coord_r[1][1];
	temp[1][0] = ch_coord_r[2][0];
	temp[1][1] = ch_coord_r[2][1];
	bov_points_update(tempDraw, temp, 2);

	for (; i < nPoints; i++) {
		while (orientation(stack_r[top - 1], stack_r[top], sorted_coord[i], coord) <= 0) {
			if (top == 1) {
				stack_r[top] = sorted_coord[i];
				ch_coord_r[top][0] = coord[sorted_coord[i]][0];
				ch_coord_r[top][1] = coord[sorted_coord[i]][1];
				bov_points_update(chDraw_r, ch_coord_r, top);
				i++;
			}
			else {
				stack_r[top] = -1;
				top--;
				temp[0][0] = ch_coord_r[top][0];
				temp[0][1] = ch_coord_r[top][1];
				temp[1][0] = coord[sorted_coord[i]][0];
				temp[1][1] = coord[sorted_coord[i]][1];
				bov_points_update(tempDraw, temp, 2);

				t = clock();
				while (!bov_window_should_close(window) && difftime(clock(), t) < pause) {
					bov_line_strip_draw_with_order(window, coordDraw, order_sort, 0, nPoints);
					bov_line_strip_draw(window, chDraw_r, 0, top + 1);
					bov_line_strip_draw(window, tempDraw, 0, 2);
					bov_points_draw(window, coordDraw, 0, nPoints);
					bov_points_draw(window, current_low_point, 0, 1);
					bov_points_draw(window, current_high_point, 0, 1);
					bov_text_draw(window, text_algo);
					bov_window_update(window);
				}
			}			
		}

		top++;
		stack_r[top] = sorted_coord[i];
		ch_coord_r[top][0] = coord[sorted_coord[i]][0];
		ch_coord_r[top][1] = coord[sorted_coord[i]][1];
		bov_points_update(chDraw_r, ch_coord_r, top);

		t = clock();
		while (!bov_window_should_close(window) && difftime(clock(), t) < pause) {
			bov_line_strip_draw_with_order(window, coordDraw, order_sort, 0, nPoints);
			bov_line_strip_draw(window, chDraw_r, 0, top + 1);
			bov_line_strip_draw(window, tempDraw, 0, 2);
			bov_points_draw(window, coordDraw, 0, nPoints);
			bov_points_draw(window, current_low_point, 0, 1);
			bov_points_draw(window, current_high_point, 0, 1);
			bov_text_draw(window, text_algo);
			bov_window_update(window);
		}
	}
	int r = top + 1;
	bov_points_update(chDraw_r, ch_coord_r, r);
	/****************** RIGHT PART ******************/

	/****************** LEFT  PART ******************/
	int* stack_l = malloc(nPoints * sizeof(int));
	if (!stack_l)
		exit(0);
	bov_text_update(text_algo, (GLubyte[]) { "Then, the left part" });
	t = clock();
	while (!bov_window_should_close(window) && difftime(clock(), t) < 1000) {
		bov_line_strip_draw_with_order(window, coordDraw, order_sort, 0, nPoints);
		bov_line_strip_draw(window, chDraw_r, 0, r);
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_points_draw(window, current_low_point, 0, 1);
		bov_points_draw(window, current_high_point, 0, 1);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}
	
	for (top = 0; top < 3; top++) {
		stack_l[top] = sorted_coord[nPoints - 1 - top];
		ch_coord_l[top][0] = coord[sorted_coord[nPoints - 1 - top]][0];
		ch_coord_l[top][1] = coord[sorted_coord[nPoints - 1 - top]][1];
	}
	top--;
	i = nPoints - 4;
	while (orientation(stack_l[0], stack_l[1], stack_l[2], coord) < 0) {
		stack_l[1] = stack_l[2];
		ch_coord_l[1][0] = ch_coord_l[2][0];
		ch_coord_l[1][1] = ch_coord_l[2][1];
		stack_l[2] = sorted_coord[i];
		ch_coord_l[2][0] = coord[sorted_coord[i]][0];
		ch_coord_l[2][1] = coord[sorted_coord[i]][1];
		bov_points_update(chDraw_l, ch_coord_l, top);
		t = clock();
		while (!bov_window_should_close(window) && difftime(clock(), t) < pause) {
			bov_line_strip_draw_with_order(window, coordDraw, order_sort, 0, nPoints);
			bov_line_strip_draw(window, chDraw_r, 0, r);
			bov_line_strip_draw(window, chDraw_l, 0, nPoints);
			bov_points_draw(window, coordDraw, 0, nPoints);
			bov_points_draw(window, current_low_point, 0, 1);
			bov_points_draw(window, current_high_point, 0, 1);
			bov_text_draw(window, text_algo);
			bov_window_update(window);
		}
		i--;
	}
	
	temp[0][0] = ch_coord_l[1][0];
	temp[0][1] = ch_coord_l[1][1];
	temp[1][0] = ch_coord_l[2][0];
	temp[1][1] = ch_coord_l[2][1];
	bov_points_update(tempDraw, temp, 2);

	for (; i >= 0; i--) {
		while (orientation(stack_l[top - 1], stack_l[top], sorted_coord[i], coord) <= 0) {
			if (top == 1) {
				stack_l[top] = sorted_coord[i];
				ch_coord_l[top][0] = coord[sorted_coord[i]][0];
				ch_coord_l[top][1] = coord[sorted_coord[i]][1];
				bov_points_update(chDraw_l, ch_coord_l, top);
				i--;
			}
			else {
				stack_l[top] = -1;
				top--;

				temp[0][0] = ch_coord_l[top][0];
				temp[0][1] = ch_coord_l[top][1];
				temp[1][0] = coord[sorted_coord[i]][0];
				temp[1][1] = coord[sorted_coord[i]][1];
				bov_points_update(tempDraw, temp, 2);

				t = clock();
				while (!bov_window_should_close(window) && difftime(clock(), t) < pause) {
					bov_line_strip_draw_with_order(window, coordDraw, order_sort, 0, nPoints);
					bov_line_strip_draw(window, chDraw_r, 0, r);
					bov_line_strip_draw(window, chDraw_l, 0, top + 1);
					bov_line_strip_draw(window, tempDraw, 0, 2);
					bov_points_draw(window, coordDraw, 0, nPoints);
					bov_points_draw(window, current_low_point, 0, 1);
					bov_points_draw(window, current_high_point, 0, 1);
					bov_text_draw(window, text_algo);
					bov_window_update(window);
				}
			}
		}

		top++;
		stack_l[top] = sorted_coord[i];
		ch_coord_l[top][0] = coord[sorted_coord[i]][0];
		ch_coord_l[top][1] = coord[sorted_coord[i]][1];
		bov_points_update(chDraw_l, ch_coord_l, top);

		t = clock();
		while (!bov_window_should_close(window) && difftime(clock(), t) < pause) {
			bov_line_strip_draw_with_order(window, coordDraw, order_sort, 0, nPoints);
			bov_line_strip_draw(window, chDraw_r, 0, r);
			bov_line_strip_draw(window, chDraw_l, 0, top + 1);
			bov_line_strip_draw(window, tempDraw, 0, 2);
			bov_points_draw(window, coordDraw, 0, nPoints);
			bov_points_draw(window, current_low_point, 0, 1);
			bov_points_draw(window, current_high_point, 0, 1);
			bov_text_draw(window, text_algo);
			bov_window_update(window);
		}
	}
	bov_points_update(chDraw_l, ch_coord_l, top + 1);
	/****************** LEFT  PART ******************/

	bov_points_set_color(chDraw_r, (GLfloat[4]) { 0.0, 1.0, 0.0, 1.0 });
	bov_points_set_color(chDraw_l, (GLfloat[4]) { 0.0, 1.0, 0.0, 1.0 });
	bov_text_update(text_algo, (GLubyte[]) { "Convex Hull found !" });
	printf("----- DONE -----");
	while (!bov_window_should_close(window)) {
		bov_line_strip_draw(window, chDraw_r, 0, nPoints);
		bov_line_strip_draw(window, chDraw_l, 0, nPoints);
		bov_points_draw(window, coordDraw, 0, nPoints);
		bov_text_draw(window, text_algo);
		bov_window_update(window);
	}

	free(sorted_coord);
	bov_text_delete(text_algo);
	bov_points_delete(current_low_point);
	bov_points_delete(current_high_point);
	bov_points_delete(current_node_point);
	bov_points_delete(chDraw_l);
	bov_points_delete(chDraw_r);
	bov_points_delete(tempDraw);
	bov_order_delete(order_sort);
	free(ch_coord_l);
	free(ch_coord_r);
	free(temp);
	free(stack_r);
	free(stack_l);

	return;
}

static void key_callback(GLFWwindow* self, int key, int scancode, int action, int mods) {
	static unsigned screenshot_nbr = 0;
	char screenshot_name[64] = "screenshot";
	bov_window_t* window = (bov_window_t*)glfwGetWindowUserPointer(self);
	if (action == GLFW_PRESS || action == GLFW_REPEAT) {
		switch (key) {
		case GLFW_KEY_ESCAPE:
			glfwSetWindowShouldClose(self, GL_TRUE);
			break;
		case GLFW_KEY_SPACE:
			if (window->running == 0) {
				window->running = 1;
			}
			else {
				window->running = 0;
			}
			break;
		case GLFW_KEY_H:
		case GLFW_KEY_K:
			if (window->help_needed == 0) {
				window->help_needed = 1;
			}
			else {
				window->help_needed = 0;
			}
			break;
		case GLFW_KEY_P:
			snprintf(screenshot_name + 10, 54, "%u.ppm", screenshot_nbr++);
			bov_window_screenshot(window, screenshot_name);
			break;
		case GLFW_KEY_UP:
			window->counter++;
			break;
		case GLFW_KEY_DOWN:
			window->counter--;
			break;
		}

		if (!points_recorded) {
			switch (key) {
			case GLFW_KEY_R:
				points_recorded = 1;
				break;
			}
		}

		else {
			switch (key) {
			case GLFW_KEY_ENTER:
				printf("Let's go !\n");
				ready = 1;
				break;
			}
		}
	}
	if (key == GLFW_KEY_ESCAPE)
		glfwSetWindowShouldClose(self, GL_TRUE);
}

static void mouse_button_callback(GLFWwindow* self, int button, int action, int mod) {
	bov_window_t* window = (bov_window_t*)glfwGetWindowUserPointer(self);

	if (button == GLFW_MOUSE_BUTTON_LEFT) {
		if (action == GLFW_PRESS) {
			glfwSetCursor(self, window->leftClickCursor);
			window->clickTime[0] = window->wtime;
			if (nPoints == nPoints_max) {
				points_recorded = 1;
				printf("The maximum number of points is reached\n");
			}
			if (!points_recorded) {
				coord[nPoints][0] = (double)window->cursorPos[0];
				coord[nPoints][1] = window->cursorPos[1];
				nPoints_supp += 1;
				printf("Point added !\n");
			}
			else {
				points_recorded = 1;
				printf("You can't add points anymore !\n");
			}
		}
		else {
			glfwSetCursor(self, NULL);
			window->clickTime[0] = -window->wtime;
		}
	}
	else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
		if (action == GLFW_PRESS) {
			glfwSetInputMode(self, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
			window->clickTime[1] = window->wtime;
		}
		else {
			glfwSetInputMode(self, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
			window->clickTime[1] = -window->wtime;
		}
	}
}

int main() {
	// initialization of parameters used in orientation()
	my_exactinit();

	// initialization of the window
	bov_window_t* window = bov_window_new(800, 800, "My ConvexHull Algo");
	bov_window_set_color(window, (GLfloat[]) { 0.9f, 0.85f, 0.8f, 1.0f });

	// give a bit of entropy for the seed of rand()
	// or it will always be the same sequence
	int seed = (int)time(NULL);
	srand(seed);
	// we print the seed so you can get the distribution of points back
	printf("seed=%d\n\n", seed);
	printf("-----------------------------------------------------\n\n");

	// welcome
	printf("HELLO\n");
	printf("  We decided to implement algorithms which find the Convex Hull of a set of points.\n");
	printf("  You just have to follow the instructions on the terminal and then on the window.\n");
	printf("-----------------------------------------------------\n");

	// choice if you want animation or not
	int animation, s;
	printf("ANIMATION OR NOT ?\n");
	printf("  - Press 1 if you want some animations (max 1000 points)\n");
	printf("  - Press 2 if you only want the final result (max 1 million points)\n");
	printf("     --> Your choice : ");
	s = scanf("%d", &animation);
	printf("-----------------------------------------------------\n");

	// points generation
	int nPoints_scan;
	printf("HOW MANY POINTS ?\n");
	switch (animation) {
	case 1:
		printf("  Choose the number of points [10 - %d] : ", nPoints_max);
		break;
	case 2:
		printf("  Choose the number of points [10 - 1000000] : ");
		break;
	}
	s = scanf("%d", &nPoints_scan);
	printf("-----------------------------------------------------\n");

	// choice of the algorithm
	int algo;
	printf("WHICH METHOD ?\n");
	printf("  Then, you have the choice between 3 algorithms to find the Convex Hull of a set of points:\n");
	printf("    - Press 1 for Jarvis March / Gift wrapping\n");
	printf("    - Press 2 for Graham Scan with points sorted according to their polar (theta) coordinates\n");
	printf("    - Press 3 for Graham Scan with points sorted according to their cartesian (y) coordinates\n");
	printf("     --> Your choice : ");
	s = scanf("%d", &algo);
	printf("-----------------------------------------------------\n");
	
	int* ch = NULL;
	switch (animation) {
	case 1:
		ready = 0;
		points_recorded = 0;
		nPoints_supp = 0;

		// points drawing
		coord = malloc(sizeof(coord[0]) * nPoints_max);
		random_points(coord, (GLsizei)nPoints_scan);

		glfwSetKeyCallback(window->self, key_callback);
		glfwSetMouseButtonCallback(window->self, mouse_button_callback);

		bov_text_t* text_click = bov_text_new((GLubyte[]) { "Click to add points\nPress R when you are ready !" }, GL_STATIC_DRAW);
		bov_text_set_fontsize(text_click, 0.1);
		bov_text_set_boldness(text_click, 0.35);
		bov_text_set_color(text_click, (GLfloat[4]) { 0.9, 0.5, 0, 1.0 });
		bov_text_set_outline_width(text_click, 0.8);
		bov_text_set_outline_color(text_click, (GLfloat[4]) { 0.0, 0.0, 0.0, 1.0 });
		bov_text_set_pos(text_click, (GLfloat[2]) { -0.9, 0.9 });

		bov_text_t* text_search = bov_text_new((GLubyte[]) { "All the points are here\nPress ENTER" }, GL_STATIC_DRAW);
		bov_text_set_fontsize(text_search, 0.1);
		bov_text_set_boldness(text_search, 0.35);
		bov_text_set_color(text_search, (GLfloat[4]) { 0.2, 0.1, 0.8, 1.0 });
		bov_text_set_outline_width(text_search, 0.8);
		bov_text_set_outline_color(text_search, (GLfloat[4]) { 0.0, 0.0, 0.0, 1.0 });
		bov_text_set_pos(text_search, (GLfloat[2]) { -0.9, 0.9 });

		bov_text_t* text_zoom = bov_text_new((GLubyte[]) { "You might have to dezoom a bit for the next step" }, GL_STATIC_DRAW);
		bov_text_set_fontsize(text_zoom, 0.08);
		bov_text_set_boldness(text_zoom, 0.35);
		bov_text_set_color(text_zoom, (GLfloat[4]) { 0.6, 0.45, 0.3, 1.0 });
		bov_text_set_outline_width(text_zoom, 0.8);
		bov_text_set_outline_color(text_zoom, (GLfloat[4]) { 0.0, 0.0, 0.0, 1.0 });
		bov_text_set_pos(text_zoom, (GLfloat[2]) { -0.96, 0 });

		// the user can add points by clicking in the window
		while (!bov_window_should_close(window) && !ready) {
			if (!points_recorded)
				bov_text_draw(window, text_click);
			else {
				bov_text_draw(window, text_search);
				bov_text_draw(window, text_zoom);
			}
			nPoints = nPoints_scan + nPoints_supp;
			coordDraw = bov_points_new(coord, nPoints, GL_STATIC_DRAW);
			bov_points_set_color(coordDraw, (GLfloat[4]) { 0.0, 0.0, 0.0, 1.0 });
			bov_points_set_width(coordDraw, 0.005);
			bov_points_set_outline_width(coordDraw, -1.0);
			bov_points_draw(window, coordDraw, 0, nPoints);
			bov_window_update(window);
		}

		if (ready) {
			bov_points_draw(window, coordDraw, 0, nPoints);
			bov_window_update(window);
		}

		// we can run the algorithm
		ch = malloc(nPoints * sizeof(int));
		if (!ch)
			exit(0);
		switch (algo) {
		case 1:
			convex_hull_jarvis_march_animation(coord, nPoints, ch, window);
			break;
		case 2:
			convex_hull_graham_scan_angle_animation(coord, nPoints, ch, window);
			break;
		case 3:
			convex_hull_graham_scan_y_animation(coord, nPoints, ch, window);
			break;
		}

		bov_text_delete(text_click);
		bov_text_delete(text_search);
		free(coord); free(coordDraw);
		free(ch);
		break;

	case 2:
		// points generation
		nPoints = (GLsizei)nPoints_scan;
		printf("generation of %d random points ... ", (int)nPoints);
		coord = malloc(sizeof(coord[0]) * nPoints);
		random_points(coord, nPoints);
		printf("done !\n");

		// points drawing
		coordDraw = bov_points_new(coord, nPoints, GL_STATIC_DRAW);
		bov_points_set_color(coordDraw, (GLfloat[4]) { 0.0, 0.0, 0.0, 1.0 });

		ch = malloc(nPoints * sizeof(int));
		if (!ch)
			exit(0);

		// find the Convex Hull
		clock_t t1 = clock();
		switch (algo) {
		case 1:
			convex_hull_jarvis_march(coord, nPoints, ch);
			break;
		case 2:
			convex_hull_graham_scan_angle(coord, nPoints, ch);
			break;
		case 3:
			convex_hull_graham_scan_y(coord, nPoints, ch);
			break;
		}
		clock_t t2 = clock();
		printf(" --> EXECUTION TIME: %1.6f [s]\n", (float)(t2 - t1) / CLOCKS_PER_SEC);

		// Place the values in a vector with the right size
		// to create the "order" object for the drawing
		int n_ch = goodSize(ch, nPoints);
		int* ch_goodSize = malloc(n_ch * sizeof(int));
		if (!ch_goodSize)
			exit(0);
		fill(ch, n_ch, ch_goodSize);
		bov_order_t* order_ch = bov_order_new(ch_goodSize, n_ch, GL_STATIC_DRAW);

		while (!bov_window_should_close(window)) {
			// display of the points
			bov_points_set_width(coordDraw, 0.005);
			bov_points_set_outline_width(coordDraw, -1.0);
			bov_points_draw(window, coordDraw, 0, nPoints);

			// display of the convex hull
			bov_fast_line_strip_draw_with_order(window, coordDraw, order_ch, 0, nPoints);

			bov_window_update(window);
		}

		free(coord);
		free(ch);
		free(ch_goodSize);
		bov_points_delete(coordDraw);
		bov_order_delete(order_ch);
		break;
	}

	bov_window_delete(window);

	return EXIT_SUCCESS;
}