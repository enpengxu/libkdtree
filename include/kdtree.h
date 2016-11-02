/* 
 * This file is part of the libkdtree distribution (https://github.com/enpengxu/libkdtree).
 * Copyright (c) 2016 Leo Xu.
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef __KDTREE_H
#define __KDTREE_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef KD_DIM
   #error "Please define KD_DIM"
#endif

#if defined KD_FLOAT
   typedef float kdfloat;
   #define kd_sqrt(x) sqrtf(x)
   #define kd_fabs(x) fabsf(x)
   #define KD_EPSILON FLT_EPSILON
#elif defined KD_DOUBLE
   typedef double kdfloat;
   #define kd_sqrt(x) sqrt(x)
   #define kd_fabs(x) fabs(x)
   #define KD_EPSILON DBL_EPSILON
#elif defined(KD_INT)
   typedef int kdfloat;
   #define kd_sqrt(x) sqrtf((float)x)
   #define kd_fabs(x) abs(x)
   #define KD_EPSILON 0
#else
   #error "Please define type, KD_FLOAT, KD_DOUBLE or KD_INT "
#endif

#include <math.h>
#include <stdint.h>

struct kdpoint {
	kdfloat val[KD_DIM];
};
struct kdpoints {
	int32_t num;
	struct kdpoint * points;
};
struct kdplane {
	int32_t axis;
	int32_t point;
};
struct kdnode {
	struct kdplane  plane;
	struct kdnode * left;
	struct kdnode * right;
};
struct kdtree {
	int32_t dim;
	int32_t depth;
	struct kdnode * root;
	struct kdpoints points;
};
struct kdnearest {
	int32_t point;
	kdfloat distance;
};

enum KD_ALGO_CTL {
	KD_ALGO_CTL_AXIS=0,
	KD_ALGO_CTL_POS,
	KD_ALGO_CTL_LAST,
};
enum KD_ALGO_LIST {
	KD_ALGO_AXIS_VARIANCE=0,
	KD_ALGO_AXIS_INTURN,
	KD_ALGO_POS_MID,
	KD_ALGO_LAST
};

#define kdnode_point(tree, node) &(tree)->points.points[((node)->plane.point)]
#define kdnode_point_val(tree, node) \
	(tree)->points.points[((node)->plane.point)].val[(node)->plane.axis]

void kdtree_init(struct kdtree * tree);
int  kdtree_build(struct kdtree * tree, struct kdpoints * points);
int  kdtree_save(struct kdtree * tree, const char * path);
int  kdtree_load(struct kdtree * tree, const char * path);
void kdtree_free(struct kdtree * tree);
int  kdtree_nearest(struct kdtree * tree,
		const struct kdpoint * point,
		struct kdnearest * result);
int kdtree_option(enum KD_ALGO_CTL ctl, enum KD_ALGO_LIST set);

#ifdef __cplusplus
}
#endif

#endif
