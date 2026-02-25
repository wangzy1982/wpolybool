#ifndef _WPOLY_
#define _WPOLY_

#if defined(_WINDOWS)
#if defined(WPOLY_EXPORTS)
#define WPOLY_API __declspec(dllexport)
#elif defined(WPOLY_STATIC)
#define WPOLY_API
#else
#define WPOLY_API __declspec(dllimport)
#endif
#else
#define WPOLY_API
#endif

#define WPOLY_API_C extern "C" WPOLY_API

/**
 * @brief Creates a new empty polygon object
 *
 * Allocates memory for a new polygon structure that can contain multiple loops.
 * The polygon must be freed with free_polygon() to avoid memory leaks.
 *
 * @return Pointer to the newly created polygon object, or NULL if memory allocation fails
 */
WPOLY_API_C void* new_polygon();

/**
 * @brief Frees the memory allocated for a polygon object
 *
 * Releases all memory associated with the polygon, including all contained loops
 * and their edges. Does nothing if the polygon pointer is NULL.
 *
 * @param polygon Pointer to the polygon object to free (can be NULL)
 */
WPOLY_API_C void free_polygon(void* polygon);

/**
 * @brief Creates a new empty loop object
 *
 * Allocates memory for a new loop structure that can contain multiple edges (lines, arcs, NURBS).
 * Loops are typically added to polygons using add_loop().
 * The loop must be freed as part of the polygon via free_polygon(), or individually if not added.
 *
 * @return Pointer to the newly created loop object, or NULL if memory allocation fails
 */
WPOLY_API_C void* new_loop();

/**
 * @brief Adds a loop to a polygon
 *
 * Attaches a loop to the specified polygon. The polygon takes ownership of the loop,
 * so the loop should not be freed individually after being added.
 *
 * @param polygon Pointer to the target polygon (must be non-NULL)
 * @param loop Pointer to the loop to add (must be non-NULL)
 * @return index of loop
 * @warning add_edge must be called before add_loop.
 */
WPOLY_API_C void add_loop(void* polygon, void* loop);

/**
 * @brief Creates a new 2D line segment edge
 *
 * Allocates memory for a line edge defined by start and end points in 2D space.
 * The edge can be added to a loop using add_edge().
 *
 * @param start_x X-coordinate of the line start point
 * @param start_y Y-coordinate of the line start point
 * @param end_x X-coordinate of the line end point
 * @param end_y Y-coordinate of the line end point
 */
WPOLY_API_C void* new_line2d_edge(double start_x, double start_y, double end_x, double end_y);

/**
 * @brief Creates a new 2D circular arc edge
 *
 * Allocates memory for an arc edge defined by center point, radius, start angle,
 * and angular extent (delta angle) in 2D space. Angles are in radians.
 * Positive delta_angle indicates counter-clockwise direction, negative indicates clockwise.
 *
 * @param center_x X-coordinate of the arc center point
 * @param center_y Y-coordinate of the arc center point
 * @param radius Radius of the arc (must be positive)
 * @param start_angle Start angle of the arc in radians (0 = positive X-axis)
 * @param delta_angle Angular extent of the arc in radians (can be positive/negative)
 * @return Pointer to the newly created arc edge, or NULL if memory allocation fails
 */
WPOLY_API_C void* new_arc2d_edge(double center_x, double center_y, double radius, double start_angle, double delta_angle);

/**
 * @brief Creates a new 2D NURBS curve edge
 *
 * Allocates memory for a NURBS (Non-Uniform Rational B-Spline) curve edge in 2D space.
 * Control points are provided as a flat array (x0, y0, x1, y1, ..., xn, yn).
 * Weights array must have the same length as the number of control points or NULL.
 *
 * @param degree Degree of the NURBS curve (must be non-negative, typically 1-15)
 * @param control_point_count Number of control points (must be > degree)
 * @param knots Array of knot values
 * @param control_points Flat array of control point coordinates (2 * control_point_count elements)
 * @param weights Array of weight values for each control point (control_point_count elements or NULL)
 * @return Pointer to the newly created NURBS edge, or NULL if memory allocation fails
 */
WPOLY_API_C void* new_nurbs2d_edge(int degree, int control_point_count, const double* knots, const double* control_points, const double* weights);

/**
 * @brief Adds an edge to a loop
 *
 * Attaches an edge (line, arc, or NURBS) to the specified loop. The loop takes ownership of the edge,
 * so the edge should not be freed individually after being added.
 *
 * @param loop Pointer to the target loop (must be non-NULL)
 * @param edge Pointer to the edge to add (must be non-NULL)
 */
WPOLY_API_C void add_edge(void* loop, void* edge);

/**
 * @brief Gets the number of loops in a polygon
 *
 * Returns the count of loops (outer contours and holes) contained in the polygon.
 *
 * @param polygon Pointer to the polygon (must be non-NULL)
 * @return Number of loops (>= 0)
 */
WPOLY_API_C int get_loop_count(void* polygon);

/**
 * @brief Retrieves a specific loop from a polygon by index
 *
 * Gets the loop at the specified zero-based index from the polygon.
 * The returned loop pointer is owned by the polygon and should not be freed individually.
 *
 * @param polygon Pointer to the polygon (must be non-NULL)
 * @param loop_index Zero-based index of the loop to retrieve (0 <= index < get_loop_count())
 * @return Pointer to the loop, or NULL if index is out of bounds or polygon is invalid
 */
WPOLY_API_C void* get_loop(void* polygon, int loop_index);

/**
 * @brief Gets the number of edges in a loop
 *
 * Returns the count of edges (lines, arcs, NURBS) contained in the loop.
 *
 * @param loop Pointer to the loop (must be non-NULL)
 * @return Number of edges (>= 0)
 */
WPOLY_API_C int get_edge_count(void* loop);

/**
 * @brief Retrieves a specific edge from a loop by index
 *
 * Gets the edge at the specified zero-based index from the loop.
 * The returned edge pointer is owned by the loop and should not be freed individually.
 *
 * @param loop Pointer to the loop (must be non-NULL)
 * @param edge_index Zero-based index of the edge to retrieve (0 <= index < get_edge_count())
 * @return Pointer to the edge, or NULL if index is out of bounds or loop is invalid
 */
WPOLY_API_C void* get_edge(void* loop, int edge_index);

/**
 * @brief Checks if an edge is a 2D line segment
 *
 * Determines whether the specified edge is of line segment type.
 *
 * @param edge Pointer to the edge to check (must be non-NULL)
 * @return true if the edge is a line segment, false otherwise (including invalid edge pointers)
 */
WPOLY_API_C bool is_line2d_edge(void* edge);

/**
 * @brief Checks if an edge is a 2D circular arc
 *
 * Determines whether the specified edge is of circular arc type.
 *
 * @param edge Pointer to the edge to check (must be non-NULL)
 * @return true if the edge is a circular arc, false otherwise (including invalid edge pointers)
 */
WPOLY_API_C bool is_arc2d_edge(void* edge);

/**
 * @brief Checks if an edge is a 2D NURBS curve
 *
 * Determines whether the specified edge is of NURBS curve type.
 *
 * @param edge Pointer to the edge to check (must be non-NULL)
 * @return true if the edge is a NURBS curve, false otherwise (including invalid edge pointers)
 */
WPOLY_API_C bool is_nurbs2d_edge(void* edge);

/**
 * @brief Retrieves the data of a 2D line segment edge
 *
 * Extracts the start and end point coordinates from a line edge.
 * Call is_line2d_edge() first to verify the edge type before calling this function.
 *
 * @param edge Pointer to the edge (must be non-NULL)
 * @param[out] start_x Output parameter for X-coordinate of the line start point
 * @param[out] start_y Output parameter for Y-coordinate of the line start point
 * @param[out] end_x Output parameter for X-coordinate of the line end point
 * @param[out] end_y Output parameter for Y-coordinate of the line end point
 */
WPOLY_API_C void get_line2d_edge_data(void* edge, double& start_x, double& start_y, double& end_x, double& end_y);

/**
 * @brief Retrieves the data of a 2D circular arc edge
 *
 * Extracts the geometric parameters from an arc edge. Angles are in radians.
 * Call is_arc2d_edge() first to verify the edge type before calling this function.
 *
 * @param edge Pointer to the edge (must be non-NULL)
 * @param[out] center_x Output parameter for X-coordinate of the arc center
 * @param[out] center_y Output parameter for Y-coordinate of the arc center
 * @param[out] radius Output parameter for the arc radius
 * @param[out] start_angle Output parameter for the arc start angle (radians)
 * @param[out] delta_angle Output parameter for the arc angular extent (radians)
 */
WPOLY_API_C void get_arc2d_edge_data(void* edge, double& center_x, double& center_y, double& radius, double& start_angle, double& delta_angle);

/**
 * @brief Retrieves the data of a 2D NURBS curve edge
 *
 * Extracts the NURBS curve parameters (degree, control points, weights).
 * The returned pointers point to internal data - do not free them.
 * Call is_nurbs2d_edge() first to verify the edge type before calling this function.
 *
 * @param edge Pointer to the edge (must be non-NULL)
 * @param[out] degree Output parameter for the NURBS curve degree
 * @param[out] control_point_count Output parameter for the number of control points
 * @param[out] knots Output parameter for pointer to knot array
 * @param[out] control_points Output parameter for pointer to flat control point array (x0,y0,x1,y1,...)
 * @param[out] weights Output parameter for pointer to weight array
 * @warning Do not modify or free the returned knots and control_points and weights pointers
 */
WPOLY_API_C void get_nurbs2d_edge_data(void* edge, int& degree, int& control_point_count, const double*& knots, const double*& control_points, const double*& weights);

#endif