#ifndef _WBOOL_
#define _WBOOL_

#include "wpoly.h"

/**
 * @brief  Compute the intersection of two polygons (Boolean AND operation)
 *
 * This function performs the intersection operation on two input polygons and returns a new polygon
 * that contains all regions belonging to both polygon0 and polygon1 simultaneously. During the operation,
 * the specified distance tolerance is used to handle floating-point precision issues, avoiding calculation
 * errors caused by minor numerical inaccuracies.
 *
 * @param  polygon0        Pointer to the first input polygon (non-NULL), the specific type is defined by the polygon library
 * @param  polygon1        Pointer to the second input polygon (non-NULL), the specific type is defined by the polygon library
 * @param  distance_epsilon Distance tolerance (precision threshold) used to determine if points/edges are coincident
 * @return void*           Returns a pointer to the new intersection polygon on success, or NULL on failure
 *                         (e.g., invalid input, no intersection area, etc.)
 * @note   1. The caller is responsible for freeing the memory of the returned polygon to avoid memory leaks
 *         2. Both input polygons must be valid and non-self-intersecting
 *         3. An overly small distance_epsilon may cause floating-point precision issues, while an overly large value
 *            may lead to unintended merging of edges/points
 */
WPOLY_API_C void* intersect(void* polygon0, void* polygon1, double distance_epsilon);

/**
 * @brief  Compute the union of two polygons (Boolean OR operation)
 *
 * This function performs the union operation on two input polygons and returns a new polygon
 * that contains all regions belonging to polygon0, polygon1, or both. During the operation,
 * the specified distance tolerance is used to handle floating-point precision issues, ensuring
 * correct merging of adjacent edges/points.
 *
 * @param  polygon0        Pointer to the first input polygon (non-NULL), the specific type is defined by the polygon library
 * @param  polygon1        Pointer to the second input polygon (non-NULL), the specific type is defined by the polygon library
 * @param  distance_epsilon Distance tolerance (precision threshold) used to determine if points/edges are coincident
 * @return void*           Returns a pointer to the new union polygon on success, or NULL on failure (e.g., invalid input)
 * @note   1. The caller is responsible for freeing the memory of the returned polygon to avoid memory leaks
 *         2. Both input polygons must be valid and non-self-intersecting
 *         3. An overly small distance_epsilon may cause floating-point precision issues, while an overly large value
 *            may lead to unintended merging of edges/points
 */
WPOLY_API_C void* unite(void* polygon0, void* polygon1, double distance_epsilon);

/**
 * @brief  Compute the difference of two polygons (Boolean subtraction operation)
 *
 * This function performs the subtraction operation on two input polygons and returns a new polygon
 * that contains all regions belonging to polygon0 but not to polygon1. During the operation,
 * the specified distance tolerance is used to handle floating-point precision issues, avoiding
 * calculation errors caused by minor numerical inaccuracies.
 *
 * @param  polygon0        Pointer to the minuend polygon (non-NULL), the specific type is defined by the polygon library
 * @param  polygon1        Pointer to the subtrahend polygon (non-NULL), the specific type is defined by the polygon library
 * @param  distance_epsilon Distance tolerance (precision threshold) used to determine if points/edges are coincident
 * @return void*           Returns a pointer to the new difference polygon on success, or NULL on failure
 *                         (e.g., invalid input, no remaining area after subtraction, etc.)
 * @note   1. The caller is responsible for freeing the memory of the returned polygon to avoid memory leaks
 *         2. Both input polygons must be valid and non-self-intersecting
 *         3. An overly small distance_epsilon may cause floating-point precision issues, while an overly large value
 *            may lead to unintended clipping of edges/points
 *         4. Subtraction operation is order-dependent: subtract(A,B) ≠ subtract(B,A)
 */
WPOLY_API_C void* subtract(void* polygon0, void* polygon1, double distance_epsilon);

#endif