#pragma once
#ifndef _MOVE_H
#define _MOVE_H

#include "TRIANGLE.h"



void trans(TRIANGLE *tr, Vec dir, double Value) {
	for (int i = 0; i < 3; i++) {
		tr->v[i].x += dir.x*Value;
		tr->v[i].y += dir.y*Value;
		tr->v[i].z += dir.z*Value;
		tr->bbox[0][0] = std::min(std::min(tr->v[0].x, tr->v[1].x), tr->v[2].x);
		tr->bbox[0][1] = std::min(std::min(tr->v[0].y, tr->v[1].y), tr->v[2].y);
		tr->bbox[0][2] = std::min(std::min(tr->v[0].z, tr->v[1].z), tr->v[2].z);
		tr->bbox[1][0] = std::max(std::max(tr->v[0].x, tr->v[1].x), tr->v[2].x);
		tr->bbox[1][1] = std::max(std::max(tr->v[0].y, tr->v[1].y), tr->v[2].y);
		tr->bbox[1][2] = std::max(std::max(tr->v[0].z, tr->v[1].z), tr->v[2].z);
		//tr->normal =(tr->n[0] + tr->n[1] + tr->n[2]) / 3; 
		tr->normal = Normalize(Cross((tr->v[1] - tr->v[0]), (tr->v[2] - tr->v[0])));
	}
}


void change_mat(TRIANGLE *tr, Material newmat) {
	tr->mat = newmat;
}

#endif _MOVE_H