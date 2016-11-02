/*
 * This file is part of the libkdtree distribution
 * (https://github.com/enpengxu/libkdtree).
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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <string.h>
#include "kdtree.h"

struct kdsel {
	int   at;
	int   num;
	int   mid; /* new selected */
};
struct kdbuild_ctx {
	int dim;
	struct kdpoints * points;
	struct kdsel sel;
	int * point_set;
	struct kdnode * parent;
};
struct kdsearch_ctx {
	int num;
	struct kdnode ** path;
	struct kdtree  * tree;
	const struct kdpoint * target;
	struct kdnearest nearest;
};
struct kdbuild_algo {
	int (*select_axis)(struct kdbuild_ctx * ctx);
	int (*select_pos)(struct kdbuild_ctx * ctx, int axis);
};

static kdfloat
kdpoint_distance(const struct kdpoint *a,
				 const struct kdpoint *b, int dim);
static void
kdpoints_free(struct kdpoints * pts);
static kdfloat
kdtree_variance(struct kdbuild_ctx * ctx,int axis);
static int
kdtree_select_axis(struct kdbuild_ctx * ctx);
static int
kdtree_select_axis_variance(struct kdbuild_ctx * ctx);
static int
kdtree_select_axis_inturn(struct kdbuild_ctx * ctx);
static int
kdtree_select_pos(struct kdbuild_ctx * ctx, int axis);
static int
kdtree_select_pos_mid(struct kdbuild_ctx * ctx, int axis);
static void
kdtree_select_hyperplane(struct kdbuild_ctx * ctx,
		struct kdplane * plane);
static void
kdtree_split(const struct kdbuild_ctx * ctx,
			 const struct kdplane * plane,
			 struct kdsel * lsel,
			 struct kdsel * rsel);
static void
kdtree_node_build(struct kdbuild_ctx * ctx,
				  struct kdnode * node);
static int
kdtree_node_depth (struct kdnode * node);
static void
kdnode_nearest(struct kdsearch_ctx * ctx);
static void
kdnode_free(struct kdnode * node);

static struct kdbuild_algo kdbuild_algo = {
	.select_axis = kdtree_select_axis_variance,
	.select_pos = kdtree_select_pos_mid,
};
int
kdtree_option(enum KD_ALGO_CTL ctl, enum KD_ALGO_LIST set)
{
	switch (ctl) {
	case KD_ALGO_CTL_AXIS:
		switch (set) {
		case KD_ALGO_AXIS_VARIANCE:
			kdbuild_algo.select_axis = kdtree_select_axis_variance;
			break;
		case KD_ALGO_AXIS_INTURN:
			kdbuild_algo.select_axis = kdtree_select_axis_inturn;
			break;
		default:
			return -EINVAL;
		}
	case KD_ALGO_CTL_POS:
		switch (set) {
		case KD_ALGO_POS_MID:
			kdbuild_algo.select_pos = kdtree_select_pos_mid;
			break;
		default:
			return -EINVAL;
	}
	default:
		return -EINVAL;
	}
}


int
kdtree_build(struct kdtree * tree,
			 struct kdpoints * points)
{
	if (!tree || !points){
		return -EINVAL;
	}
	if (!points->points || points->num<=0) {
		return -EINVAL;
	}
	tree->dim = KD_DIM;
	tree->depth = 0;
	tree->root = 0;
	tree->points = *points;

	int i, * set;
	set = (int *)malloc(sizeof(int)*points->num);
	if (!set)
		return -ENOMEM;

	for (i=0; i<points->num; i++) {
		set[i] = i;
	}

	tree->root = (struct kdnode *)malloc(sizeof(*tree->root));
	if (!tree->root){
		free(set);
		return -ENOMEM;
	}
	tree->root->left = 0;
	tree->root->right= 0;

	struct kdbuild_ctx ctx = {
		.dim = tree->dim,
		.points = points,
		.sel = {
			.at = 0,
			.num = points->num
		},
		.point_set = set,
		.parent = 0,
	};
	kdtree_node_build(&ctx, tree->root);
	free(set);
	tree->depth = kdtree_node_depth(tree->root);
	assert(tree->depth>0);
	return 0;
}
void
kdtree_free(struct kdtree * tree)
{
	if (!tree)
		return ;
	kdnode_free(tree->root);
	kdpoints_free(&tree->points);
	tree->root = 0;
}
void
kdtree_init(struct kdtree * tree)
{
	tree->root = 0;
	tree->points.num = 0;
	tree->points.points = 0;
	tree->depth = 0;
	tree->dim = KD_DIM;
}

#define kdsearch_add_node(ctx, n) \
	(ctx)->path[(ctx)->num++] = n; \
	assert((ctx)->num <=(ctx)->tree->depth*2)

int
kdtree_nearest(struct kdtree * tree,
			   const struct kdpoint * point,
			   struct kdnearest * result)
{
	struct kdnode * node = tree->root;
	struct kdsearch_ctx ctx;
	ctx.num = 0;
	ctx.target = point;
	ctx.tree = tree;
	ctx.nearest.point = -1;

	assert(tree->depth>0);
	/*TODO. put a cache buffer in tree to avoid alloc memory each time */
	ctx.path = (struct kdnode **)malloc(sizeof(struct kdnode *)*tree->depth*2);
	if (!ctx.path){
		return -ENOMEM;
	}
	while (node) {
		kdsearch_add_node(&ctx, node);
		kdfloat det = point->val[node->plane.axis] -
			kdnode_point_val(tree, node);
		if (det < 0)
			node = node->left;
		else
			node = node->right;
	}
	/* begin to search */
	kdnode_nearest(&ctx);

	free(ctx.path);
	*result = ctx.nearest;
	return 0;
}

static void
kdtree_node_build(struct kdbuild_ctx * ctx,
	struct kdnode * node)
{
	if (!ctx->sel.num)
		return ;

	kdtree_select_hyperplane(ctx, &node->plane);
	if (ctx->sel.num == 1) {
		return;
	}

	struct kdsel lsel, rsel;
	kdtree_split(ctx, &node->plane, &lsel, &rsel);

	if (lsel.num > 0){
		struct kdbuild_ctx lctx = *ctx;
		lctx.sel = lsel;
		lctx.parent = node;

		node->left = (struct kdnode *)malloc(sizeof(struct kdnode));
		node->left->left  = 0;
		node->left->right = 0;

		kdtree_node_build(&lctx, node->left);
	}
	if (rsel.num > 0){
		struct kdbuild_ctx rctx = *ctx;
		rctx.sel = rsel;
		rctx.parent = node;

		node->right = (struct kdnode *)malloc(sizeof(struct kdnode));
		node->right->left  = 0;
		node->right->right = 0;
		kdtree_node_build(&rctx, node->right);
	}
}

static int
kdtree_node_depth(struct kdnode * node)
{
	if (!node)
		return 0;
	int dl = 0, dr = 0;
	if (node->left)
		dl = kdtree_node_depth(node->left);
	if (node->right)
		dr = kdtree_node_depth(node->right);
	return 1 + (dl > dr ? dl : dr);
}

static int
kdtree_select_axis_inturn(struct kdbuild_ctx * ctx)
{
	if (!ctx->parent)
		return 0;
	return (ctx->parent->plane.axis+1) % ctx->dim;
}
static int
kdtree_select_axis_variance(struct kdbuild_ctx * ctx)
{
	int i, axis;
	kdfloat vmax;
	axis = 0;
	vmax = kdtree_variance(ctx, 0);
	for(i=1; i<ctx->dim; i++){
		kdfloat v = kdtree_variance(ctx, i);
		if (v > vmax){
			axis = i;
			vmax = v;
		}
	}
	return axis;
}

static int
kdtree_select_axis(struct kdbuild_ctx * ctx)
{
	struct kdsel * sel = &ctx->sel;
	assert(sel->num > 1);
	/*TODO. algorithm functinor */
	return kdbuild_algo.select_axis(ctx);
}

struct kdtree_qort_ctx {
	int axis;
	struct kdpoints * points;
};


static int
kdtree_cmp(const void *a, const void *b, void *ctx)
{
	struct kdtree_qort_ctx * c = (struct kdtree_qort_ctx *)ctx;

	int ia = *(int *)a;
	int ib = *(int *)b;

	assert(ia >=0 && ia < c->points->num);
	assert(ib >=0 && ib < c->points->num);

	kdfloat va = c->points->points[ia].val[c->axis];
	kdfloat vb = c->points->points[ib].val[c->axis];

	if (va < vb)
		return -1;
	if (va > vb)
		return 1;
	return 0;
}

static void
kdtree_sort(struct kdpoints * points,
			int * set, int start, int num, int axis)
{
	/* _axis = axis; */
	/* qsort(&set[start], num, sizeof(int), kdtree_cmp); */
	struct kdtree_qort_ctx c = {
		.axis = axis,
		.points = points,
	};
#if (defined __linux__)
	qsort_r(&set[start], num, sizeof(int), kdtree_cmp, &c);
#elif (defined _WIN32)
	qsort_s(&set[start], num, sizeof(int), kdtree_cmp, &c);
#else
	/*TODO. other platforms */
#endif
}

static int
kdtree_select_pos_mid(struct kdbuild_ctx * ctx, int axis)
{
	/* sort array & get the middle pos */
	kdtree_sort(ctx->points, ctx->point_set,
				ctx->sel.at, ctx->sel.num, axis);
	ctx->sel.mid = ctx->sel.at + ctx->sel.num/2;
	return ctx->point_set[ctx->sel.mid];
}

static int
kdtree_select_pos(struct kdbuild_ctx * ctx, int axis)
{
	struct kdsel * sel = &ctx->sel;
	assert(sel->num > 1);
	return kdbuild_algo.select_pos(ctx, axis);
}


static void
kdtree_split(const struct kdbuild_ctx * ctx,
			 const struct kdplane * plane,
			 struct kdsel * lsel,
			 struct kdsel * rsel)
{
	const struct kdsel * sel = &ctx->sel;
	assert(sel->num > 1);
	lsel->at  = sel->at;
	lsel->num = sel->mid - sel->at;

	rsel->at = sel->mid+1;
	rsel->num = sel->num - 1 - lsel->num;

	assert(lsel->num >=0 && rsel->num >=0);
}

static void
kdtree_select_hyperplane(struct kdbuild_ctx * ctx,
		struct kdplane * plane)
{
	struct kdsel * sel = &ctx->sel;
	assert(sel->num > 0);
	if (sel->num == 1){
		plane->axis = 0;

		int * set = ctx->point_set;
		assert(sel->at < ctx->points->num);
		int idx = set[sel->at];
		assert(idx < ctx->points->num);
		plane->point = idx;
		return ;
	}
	/* select algorithm to pick up plane */
	plane->axis  = kdtree_select_axis(ctx);
	plane->point = kdtree_select_pos(ctx, plane->axis);
}

static kdfloat
kdtree_variance(struct kdbuild_ctx * ctx,int axis)
{
	int i, idx;
	kdfloat s, r, m;
	int * set = ctx->point_set;
	struct kdsel * sel = &ctx->sel;
	struct kdpoint * pts = ctx->points->points;
	for (s=0, i=0; i<sel->num; i++) {
		idx = sel->at + i;
		assert(idx<ctx->points->num);
		idx = set[idx];
		assert(idx<ctx->points->num);

		s += pts[idx].val[axis];
	}
	s /= sel->num;
	for (r=0, i=0; i<sel->num; i++) {
		idx = sel->at + i;
		assert(idx<ctx->points->num);
		idx = set[idx];
		assert(idx<ctx->points->num);

		m = pts[idx].val[axis] -s;
		r += m * m;
	}
	return r;
}

static void
kdnode_free(struct kdnode * node)
{
	struct kdnode * l, *r;
	if (!node)
		return;
	l = node->left;
	r = node->right;
	free(node);
	if (l) kdnode_free(l);
	if (r) kdnode_free(r);
}
static void
kdnode_nearest(struct kdsearch_ctx * ctx)
{
	if (!ctx->num)
		return ;
	/* pop the last node */
	struct kdnode * node;
	node = ctx->path[ctx->num-1];
	ctx->num-- ;

	if (ctx->nearest.point == -1) {
		/* pickup the leaf node as nearest point */
		ctx->nearest.point = node->plane.point;
		ctx->nearest.distance = kdpoint_distance(ctx->target,
				kdnode_point(ctx->tree, node), ctx->tree->dim);
	} else {
		kdfloat dis = kdpoint_distance(ctx->target,
				kdnode_point(ctx->tree, node), ctx->tree->dim);

		if (dis < ctx->nearest.distance) {
			/* new nearest point */
			ctx->nearest.distance = dis;
			ctx->nearest.point = node->plane.point;
		}
		/* check the other side of current node */
		dis = ctx->target->val[node->plane.axis] -
			kdnode_point_val(ctx->tree, node);
		if (dis < 0){
			/* target is @ left side */
			dis = -dis;
			if (node->right && dis <= ctx->nearest.distance){
				//ctx->path[ctx->num++] = node->right;
				kdsearch_add_node(ctx, node->right);
			}
			if (node->left) {
				kdsearch_add_node(ctx, node->left);
			}
		} else {
			/* target is @ right side */
			if (node->left && dis <= ctx->nearest.distance){
				//ctx->path[ctx->num++] = node->left;
				kdsearch_add_node(ctx, node->left);
			}
			if (node->right){
				kdsearch_add_node(ctx, node->right);
			}
		}
	}
	kdnode_nearest(ctx);
}

static kdfloat
kdpoint_distance(const struct kdpoint *a,
				 const struct kdpoint *b,
				 int dim)
{
	int i;
	kdfloat r = 0;
	for (i=0; i<dim; i++) {
		r += (a->val[i] - b->val[i]) * (a->val[i] - b->val[i]);
	}
	return kd_sqrt(r);
}

static void
kdpoints_free(struct kdpoints * pts)
{
	if (pts && pts->points){
		free(pts->points);
		pts->points = 0;
		pts->num = 0;
	}
}

#define KD_MAGIC    0x40490fda
#define KD_VERSION  1
/*
 * data layout:
 * 1. header block
 * 2. node block
 * 3. point block
 * ...
 */
struct kdio_header {
	char    kdtree[4]; /* kdtr */
	int32_t magic;
	int32_t version;
	int32_t dim;
	int32_t depth;
	int32_t point_num;
	int64_t node_off;
	int64_t point_off;
	int32_t reserved[64];
};
struct kdio_node {
	struct kdplane plane;
	int32_t left; /* 0 means NULL,the first node in file is always as reserved*/
	int32_t right;
};

static int
kdnode_to_ionode(struct kdnode * node,
				 struct kdio_node * ionodes,
				 int * num, int total)
{
	if (!node)
		return 0;
	int n = *num;
	(void)total;
	assert(n>0 && n<total);

	ionodes[n].plane = node->plane;
	(*num) ++;

	ionodes[n].left  = kdnode_to_ionode(node->left,  ionodes, num, total);
	ionodes[n].right = kdnode_to_ionode(node->right, ionodes, num, total);

	assert (ionodes[n].left  >=0 && ionodes[n].left  < total);
	assert (ionodes[n].right >=0 && ionodes[n].right < total);
	return  n;
}
static struct kdnode *
ionode_to_kdnode(struct kdio_node * ionodes, int start, int total)
{
	struct kdnode * node;
	assert(start>=0 && start <total);
	if (start == 0)
		return NULL;
	node = (struct kdnode *)malloc(sizeof(*node));
	if (!node) {
		assert (0 && "no enough memory !");
		/*TODO. should exit*/
		return NULL;
	}
	node->plane = ionodes[start].plane;
	node->left  = ionode_to_kdnode(ionodes, ionodes[start].left, total);
	node->right = ionode_to_kdnode(ionodes, ionodes[start].right, total);
	return node;
}

int
kdtree_save(struct kdtree * tree, const char * path)
{
	struct kdio_header header;
	header.kdtree[0] = 'k';
	header.kdtree[1] = 'd';
	header.kdtree[2] = 't';
	header.kdtree[3] = 'r';
	header.magic = KD_MAGIC;
	header.version = KD_VERSION;
	header.dim = tree->dim;
	header.depth = tree->depth;;
	header.point_num = tree->points.num;
	header.node_off = 0;
	header.point_off = 0;
	memset(header.reserved, 0, sizeof(header.reserved));

	FILE * pf = fopen(path, "wb");
	if (!pf) {
		return -errno;
	}
	if (sizeof(header) !=
		fwrite(&header, 1, sizeof(header), pf)){
		fclose(pf);
		return -errno;
	}

	struct kdio_node * ionodes;
	size_t ionode_size = sizeof(*ionodes)*(tree->points.num+1);
	ionodes = (struct kdio_node *)malloc(ionode_size);
	if (!ionodes){
		fclose(pf);
		return -ENOMEM;
	}
	int num = 1;
	kdnode_to_ionode(tree->root, ionodes, &num, tree->points.num+1);
	header.node_off = ftell(pf);
	if (ionode_size != fwrite(ionodes, 1, ionode_size, pf)){
		fclose(pf);
		free(ionodes);
		return -errno;
	}

	size_t point_size = sizeof(struct kdpoint)*tree->points.num;
	header.point_off = ftell(pf);
	if (point_size != fwrite(tree->points.points, 1, point_size, pf)) {
		fclose(pf);
		free(ionodes);
		return -errno;
	}

	if (0 != fseek(pf, 0L, SEEK_SET)){
		fclose(pf);
		free(ionodes);
		return -errno;
	}
	if (sizeof(header) != fwrite(&header, 1, sizeof(header), pf)){
		fclose(pf);
		free(ionodes);
		return -errno;
	}
	free(ionodes);
	fclose(pf);
	return 0;
}
int
kdtree_load(struct kdtree * tree, const char * path)
{
	FILE * pf = fopen(path, "rb");
	if (!pf) {
		return -errno;
	}
	struct kdio_header header;
	if (sizeof(header) !=
		fread(&header, 1, sizeof(header), pf)){
		fclose(pf);
		return -errno;
	}
	if (header.kdtree[0] != 'k' ||
		header.kdtree[1] != 'd' ||
		header.kdtree[2] != 't' ||
		header.kdtree[3] != 'r' ||
		header.magic != KD_MAGIC ||
		header.version > KD_VERSION ||
		header.dim != KD_DIM ||
		header.point_num <= 0){
		fclose(pf);
		return -EINVAL;
	}
	struct kdio_node * ionodes;
	size_t ionode_size = sizeof(*ionodes)*(header.point_num+1);
	ionodes = (struct kdio_node *)malloc(ionode_size);
	if (!ionodes){
		fclose(pf);
		return -ENOMEM;
	}
	if (0 != fseek(pf, header.node_off, SEEK_SET)){
		fclose(pf);
		free(ionodes);
		return -errno;
	}
	if (ionode_size != fread(ionodes, 1, ionode_size, pf)){
		free(ionodes);
		fclose(pf);
		return -errno;
	}
	if (0 != fseek(pf, header.point_off, SEEK_SET)){
		fclose(pf);
		free(ionodes);
		return -errno;
	}
	size_t point_size = sizeof(struct kdpoint)*header.point_num;
	struct kdpoint * pts = (struct kdpoint *)malloc(point_size);
	if (!pts){
		free(ionodes);
		fclose(pf);
		return -ENOMEM;
	}
	if (point_size != fread(pts, 1, point_size, pf)) {
		free(pts);
		fclose(pf);
		return -errno;
	}
	kdtree_free(tree);
	tree->dim = header.dim;
	tree->depth = header.depth;
	tree->root = ionode_to_kdnode(ionodes, 1, header.point_num+1);
	tree->points.num = header.point_num;
	tree->points.points = pts;
	return 0;
}
