#ifndef _INTERSECTION_H_
#define _INTERSECTION_H_

#include "vec.h"
#include "common.h"
#include "TRIANGLE.h"

	struct Hitpoint {
		double distance;
		Vec normal;
		Vec orienting_normal;
		Vec position;
		

		Hitpoint() : distance(INF), normal(), orienting_normal(), position() {}
	};

	struct Intersection {
		Hitpoint hitpoint;
		int obj_id;
		int tri_num;
		Intersection():obj_id(-1),tri_num(-1)  {}
	};



#endif