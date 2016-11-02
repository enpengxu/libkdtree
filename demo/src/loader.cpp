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
#include <stdio.h>
#include <stdlib.h>

bool
load_bunny(const char * path, float ** pts, int *num,
		   float *xmin, float *xmax,
		   float *ymin, float *ymax,
		   float *zmin, float *zmax )
{
	int i=0, rc;
	FILE * pf = fopen(path, "r");
	if (!pf) {
		return false;
	}
	*num = 0;
	rc = fscanf(pf, "%d", num);
	if (rc <= 0 || num<=0){
		fclose(pf);
		return false;
	}
	*pts = (float *)malloc(sizeof(float)*(*num)*3);
	while(!feof(pf) && i<(*num)){
		rc = fscanf(pf, "%f %f %f",
				&(*pts)[3*i], &(*pts)[3*i+1], &(*pts)[3*i+2]);
		if(rc <= 0)
			break;
		i ++;
	}
	fclose(pf);

	if (i <(*num)){
		return false;
	}
	*xmin = (*pts)[0];
	*xmax = (*pts)[0];
	*ymin = (*pts)[1];
	*ymax = (*pts)[1];
	*zmin = (*pts)[2];
	*zmax = (*pts)[2];

#define update_bound(a1, a2, c) \
		if((c)>(a2))			\
			(a2) = (c);			\
		if((c)<(a1))			\
			(a1) = (c)

	for(i=0; i<(*num); i++){
		update_bound(*xmin, *xmax,(*pts)[3*i]);
		update_bound(*ymin, *ymax,(*pts)[3*i+1]);
		update_bound(*zmin, *zmax,(*pts)[3*i+2]);
	}
#undef update_bound
	return i == (*num) ? true: false;
}
