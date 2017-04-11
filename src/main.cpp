#define TINYOBJLOADER_IMPLEMENTATION
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <memory>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <atomic>
#include <mutex>
#include<list>
#include<vector>
#include<algorithm>
#include <iomanip>
#include "common.h"
#include "dirs.h"
#include "parallel.h"
#include "vec.h"
#include "ray.h"
#include "material.h"
#include "sphere.h"
#include "random.h"
#include "TRIANGLE.h"
#include"tiny_obj_loader.h"
#include "bvh.h"






// �}���̃��X�g
//Medium milkMedium(Color(4.5513, 5.8294, 7.136), Color(0.0015333, 0.0046, 0.019933), 1.3);
Medium milkMedium (Color (0.45513, 0.58294, 0.7136), Color(0.00015333, 0.00046, 0.0019933), 1.3);
Medium cokeMedium(Color(8.9053e-05, 8.372e-05, 0) * 1.0, Color(0.10014, 0.16503, 0.2468) * 1.0, 1.3);
Medium diamondMedium(Color(EPS,EPS,EPS) * 1.0, Color(EPS,EPS,EPS) * 1.0, 1.5);
// �}�e���A���̃��X�g
Material lightMat(Color(36.0), Color(0.0), DIFFUSE);//TODO��߂�
Material grayMat(Color(0.0), Color(0.75), DIFFUSE);
Material redMat(Color(0.0), Color(0.75, 0.25, 0.25), DIFFUSE);
Material greenMat(Color(0.0), Color(0.25, 0.75, 0.25), DIFFUSE);
Material blueMat(Color(0.0), Color(0.25, 0.25, 0.75), DIFFUSE);
Material mirrorMat(Color(0.0), Color(0.999), SPECULAR);
Material glassMat(Color(0.0), Color(0.999), REFRACTION);
Material milkMat(Color(0.0), Color(0.999), TRANSLUCENT, milkMedium);
Material cokeMat(Color(0.0), Color(0.999), TRANSLUCENT, cokeMedium);
Material diamondMat(Color(0.0), Color(0.999), REFRACTION,diamondMedium);

// �V�[���̕⏕���
const int LightID = 0;
std::vector<std::string> objList =
{
	"bunny.obj"
};
std::vector<std::vector <TRIANGLE>>triangles;

//�V�[���Ƃ̌�������֐�(�O�p�`ver)
inline bool intersect_scene_triangle(const Ray &ray, Intersection  *intersection) {
	
	intersection->hitpoint.distance = INF;

	for(int j=0;j<objList.size();j++){
		const int n = triangles[j].size();

		for (int i = 0; i < int(n); i++) {

			Hitpoint hitpoint;
			if (triangles[j][i].intersect(ray, &hitpoint)) {
				if (hitpoint.distance <= intersection->hitpoint.distance) {
					intersection->hitpoint = hitpoint;
					intersection->obj_id = j;
					intersection->tri_num = i;
					intersection->hitpoint.orienting_normal = Dot(intersection->hitpoint.normal, ray.dir) < 0.0 ? intersection->hitpoint.normal : (-1.0 * intersection->hitpoint.normal);
				}
			}
		}
	}
	return ((intersection->hitpoint.distance<INF));
}




//���̌��������(SPECULAR,REFRACTION)�̂Ƃ��g��//���Ӂ@�����̂ǂ�����inteersect�ŕK���n������:���Ă��Ȃ���Ύg��Ȃ� ref)SPECULAR

Color direct_radiance(const Vec &v0, const Vec &normal, const int id,const int tri_num, const Vec &light_pos,Intersection lintersect) {
	const Vec light_normal = triangles[lintersect.obj_id][lintersect.tri_num].normal;
	const Vec light_dir = Normalize(light_pos - v0);
	const double dist2 = (light_pos - v0).LengthSquared();
	const double dot0 = Dot(normal, light_dir);
	const double dot1 = Dot(light_normal, -1.0 * light_dir);

	if (dot0 >= 0.0 && dot1 >= 0.0) {
		const double G = dot0 * dot1 / dist2;
		Intersection intersection;
		intersect_scene_triangle(Ray(v0, light_dir), &intersection);
		if (std::abs(sqrt(dist2) - intersection.hitpoint.distance) < 1e-3) {
			Vec edge1 = (triangles[LightID][1].v[1] - triangles[LightID][1].v[0]);
			Vec edge2 = (triangles[LightID][1].v[2] - triangles[LightID][1].v[0]);
			double S = Cross(edge1,edge2).Length()/2;//light�̃��b�V����̕\�ʐ�
			return Multiply(triangles[id][tri_num].mat.ref, triangles[LightID][lintersect.tri_num].mat.Le) * (1.0 / PI) * G / (1.0 / (2*S));//TODO 2���֐��Ƃ��ċ��߂���悤�ɂ���
		
			//(1.0 / PI)�ĂȂ�ł�����
		}
	}
	return Color(0.0);
}



// ������̓_���T���v�����O���Ē��ڌ����v�Z����B//DIFFUSE�ʂŗp����
Color direct_radiance_sample(const Vec &v0, const Vec &normal, const int obj_id,const int tri_num, double u0, double u1,double u2) {
	// ������̈�_���T���v�����O����
	//u0,u1��0�`1�̗���
	int r1 = (int)(u0*(2));//TODO ���b�V���̐����������܂Ȃ��Ă��ł���悤�ɂ���.
	TRIANGLE light_triangle =triangles[LightID][r1];//�����_����light�̃��b�V������ꂽ
	Vec light_pos = light_triangle.v[0] + u1*(light_triangle.v[1] - light_triangle.v[0]) + u2*(1.0-u1)*(light_triangle.v[2] - light_triangle.v[1]);
	Vec dir = Normalize(light_pos - v0);
	Intersection lintersect;
	intersect_scene_triangle(Ray(v0, dir), &lintersect);
	if (lintersect.obj_id == LightID) {
		return direct_radiance(v0, normal, obj_id, tri_num, light_pos, lintersect);
	}
	return Color(0.0);
}

//// �w�肵���ʒu�ւ̒��ڌ����v�Z����B�������A��Ԓ��̓_�������󂯂�B
Color direct_radiance_media(const Vec &v0, const Vec &light_pos,const TRIANGLE &light_triangle) {
	const Vec light_normal = light_triangle.normal;
	const Vec light_dir = Normalize(light_pos - v0);
	const double dist2 = (light_pos - v0).LengthSquared();
	const double dot1 = Dot(light_normal, -1.0 * light_dir);

	if (dot1 >= 0) {

		Vec edge1 = (light_triangle.v[1] - light_triangle.v[0]);
		Vec edge2 = (light_triangle.v[2] - light_triangle.v[0]);
		double S = Cross(edge1, edge2).Length() / 2;//light�̃��b�V����̕\�ʐ�



		const double G = dot1 / dist2;
		Intersection intersect;
		intersect_scene_triangle(Ray(v0, light_dir), &intersect);
		if (fabs(sqrt(dist2) - intersect.hitpoint.distance) < 1e-3) {
			const Vec ret = light_triangle.mat.Le * (1.0 / PI) * G / (1.0 / 2*S);
			return ret;
		}
	}
	return Color(0.0);
}

// �}�����ɂ���_v0�ɑ΂��钼�ڌ��̉e�����v�Z����
Color direct_radiance_sample_media(const Vec &v0, const Medium &mat, double u0, double u1,double u2) {
	// ������̈�_���T���v�����O����
	int r1 = (int)(u0*(2));//TODO ���b�V���̐����������܂Ȃ��Ă��ł���悤�ɂ���.
	TRIANGLE light_triangle =triangles[LightID][r1];//�����_����light�̃��b�V������ꂽ
	Vec light_pos = light_triangle.v[0] + u1*(light_triangle.v[1] - light_triangle.v[0]) + u2*(1.0-u1)*(light_triangle.v[2] - light_triangle.v[1]);
	Vec dir = Normalize(light_pos - v0);

	Intersection lintersect;
	intersect_scene_triangle(Ray(v0, dir), &lintersect);
	
	// �}���̊O�ɏo��܂ł̌��̌������v�Z
	

	const Color sigT = mat.sigS + mat.sigA;
	const Vec transmittance_ratio = Vec::exp(-sigT * lintersect.hitpoint.distance);

	// �O�ɏo��_�ւ̒��ڌ��̉e�����v�Z����
	const Vec v1 = v0 + lintersect.hitpoint.distance * dir;
	const Color direct_light = direct_radiance_media(v1, light_pos,light_triangle);
	return Multiply(transmittance_ratio, direct_light);
}



// ray��������̕��ˋP�x�����߂�
// �{�����[�������_�����O�������Ɋ�Â�
Color radiance(const Ray &ray, const Medium &medium, Random &rng, int depth, int maxDepth) {
	// ��������
	//double t;
	//int id;
	////if (!intersect_scene(ray, &t, &id)) {
	////	return Color(0.0);
	////}



	Intersection intersection;
	if (!intersect_scene_triangle(ray, &intersection)) {
		return Color(0.0);
	}
	


	const double t = intersection.hitpoint.distance;//ok
	const int obj_id = intersection.obj_id;
	const int tri_num = intersection.tri_num;
	const TRIANGLE &obj = triangles[obj_id][tri_num];
	const Vec normal = intersection.hitpoint.normal;//ok
	const Vec orienting_normal = Dot(intersection.hitpoint.normal, ray.dir) < 0.0 ? intersection.hitpoint.normal : (-1.0 * intersection.hitpoint.normal); // �����ʒu�̖@���i���̂���̃��C�̓��o���l���j
	const Vec hitpoint = ray.org +(intersection.hitpoint.distance)* ray.dir;//ok

	return obj.mat.ref/t;

	// �ő��bounce����]��
	if (depth >= maxDepth) {
		return obj.mat.Le;
	}

	// ���ȏヌ�C��ǐՂ����烍�V�A�����[���b�g�����s���ǐՂ�ł��؂邩�ǂ����𔻒f����
	double russian_roulette_probability = std::max(obj.mat.ref.x, std::max(obj.mat.ref.y, obj.mat.ref.z));
	russian_roulette_probability = std::max(0.05, russian_roulette_probability);
	if (rng.next01() > russian_roulette_probability) {
		return obj.mat.Le / (1.0 - russian_roulette_probability);
	}

	// ���C�ƕ��̂̌����_����̕��ˋP�x���v�Z���邩�A�r���̓_�ɂ�������͂���̉e�����v�Z���邩��I�����郍�V�A�����[���b�g
	const Color sigT = medium.sigS + medium.sigA;
	const double tr_average = (sigT.x + sigT.y + sigT.z) / 3.0;
	const double sc_average = (medium.sigS.x + medium.sigS.y + medium.sigS.z) / 3.0;
	double scattering_probability = std::max(0.0, std::min(sc_average, 0.90));
	if (sigT.isZero()) {
		scattering_probability = 0.0;
	}

	if (rng.next01() < scattering_probability) {
		//���ώ��R�H��
		const double u = rng.next01();
		const double d = -log(1.0 + u * (exp(-tr_average * t) - 1.0)) / tr_average;
		const double pdf = exp(-tr_average * d) * (-tr_average / (exp(-tr_average * t) - 1.0));

		//�����C�̕��������ʏォ���l�T���v�����O
		double r1 = 2.0 * PI * rng.next01();
		double r2 = 1.0 - 2.0 * rng.next01();
		const Vec next_dir(sqrt(1.0 - r2*r2) * cos(r1), sqrt(1.0 - r2*r2) * sin(r1), r2);
		const Ray next_ray(ray.org + d * ray.dir, next_dir);

		const Vec transmittance_ratio = Vec::exp(-sigT * d);
		const Vec direct_light = direct_radiance_sample_media(next_ray.org, medium, rng.next01(), rng.next01(),rng.next01());

		// �ʑ��֐�
		const double cosTheta = Dot(-Normalize(ray.dir), Normalize(next_dir));
		const double g = 0.0;
		const double denom = std::sqrt(std::max(0.0, 1.0 + g * g - 2.0 * g * cosTheta));
		const double phase = (1.0 - g * g) / (4.0 * PI * denom * denom * denom);

		if (pdf == 0.0) {
			return Color(0.0);
		}
		else {
			const double beta = (4.0 * PI * phase) / (pdf * scattering_probability * russian_roulette_probability);
			const Color phi = Multiply(transmittance_ratio, Multiply(medium.sigS, radiance(next_ray, medium, rng, depth + 1, maxDepth)));
			return beta * phi;
		}
	}
	else {
		// ���C�ƕ��̂̌����_����̕��ˋP�x�`�B���v�Z
		const Vec transmittance_ratio = Vec::exp(-sigT * t);
		switch (obj.mat.type) {
		case DIFFUSE: {
			// ���ڌ��̃T���v�����O���s��
			if (obj_id != LightID) {
				const int shadow_ray = 1;
				Vec direct_light;
				for (int i = 0; i < shadow_ray; i++) {
					direct_light = direct_light + direct_radiance_sample(hitpoint, orienting_normal, obj_id,tri_num, rng.next01(), rng.next01(),rng.next01()) / shadow_ray;
				}
			


				// orienting_normal�̕�������Ƃ������K�������(w, u, v)�����B���̊��ɑ΂��锼�����Ŏ��̃��C���΂��B
				Vec w, u, v;
				w = orienting_normal;
				if (std::abs(w.x) > 0.1)
					u = Normalize(Cross(Vec(0.0, 1.0, 0.0), w));
				else
					u = Normalize(Cross(Vec(1.0, 0.0, 0.0), w));
				v = Cross(w, u);

				// �R�T�C�������g�����d�_�I�T���v�����O
				const double r1 = 2.0 * PI * rng.next01();
				const double r2 = rng.next01();
				const double r2s = sqrt(r2);
				Vec dir = Normalize((u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1.0 - r2)));

				// �������l���������v�Z
				return direct_light + Multiply(transmittance_ratio, Multiply(obj.mat.ref, radiance(Ray(hitpoint, dir), medium, rng, depth + 1, maxDepth))) / (1.0 - scattering_probability) / russian_roulette_probability;
			}
			else if (depth == 0) {
				return obj.mat.Le;
			}
			else {
				return Color(0.0);
			}
		} break;

		case SPECULAR: {
			
			// ���S���ʂȂ̂Ń��C�̔��˕����͌���I�B
			// ���V�A�����[���b�g�̊m���ŏ��Z����̂͏�Ɠ����B
			Intersection lintersect;//���ˌ��̏��
			Ray reflection_ray = Ray(hitpoint, ray.dir - normal * 2.0 * Dot(normal, ray.dir));
			intersect_scene_triangle(reflection_ray, &lintersect);
			Vec direct_light;
			
			if (lintersect.obj_id == LightID) {
				//direct_light = direct_radiance(hitpoint, orienting_normal, reflection_ray.org + lintersect.hitpoint.distance*reflection_ray.dir);
				direct_light = direct_radiance(hitpoint, orienting_normal, obj_id, tri_num, reflection_ray.org + lintersect.hitpoint.distance * reflection_ray.dir, lintersect);
			}
			return direct_light +
			Multiply(transmittance_ratio, radiance(Ray(hitpoint, ray.dir - normal * 2.0 * Dot(normal, ray.dir)), medium, rng, depth + 1, maxDepth)) / (1.0 - scattering_probability) / russian_roulette_probability;
		} break;
		case REFRACTION: {			Ray reflection_ray = Ray(hitpoint, ray.dir - normal * 2.0 * Dot(normal, ray.dir));

			// ���˕�������̒��ڌ��T���v�����O����
			Intersection lintersect;
			intersect_scene_triangle(reflection_ray, &lintersect);
			Vec direct_light;
			if (lintersect.obj_id == LightID) {
				direct_light = direct_radiance(hitpoint, orienting_normal, obj_id, tri_num, reflection_ray.org + lintersect.hitpoint.distance * reflection_ray.dir, lintersect);
			}
			bool into = Dot(normal, orienting_normal) > 0.0; // ���C���I�u�W�F�N�g����o��̂��A����̂��A�o��Ȃ��false

															 // Snell�̖@��
			const double nc = 1.0; // �^��̋��ܗ�
			const double nt = 1.3; // �I�u�W�F�N�g�̋��ܗ�
			const double nnt = into ? nc / nt : nt / nc;
			const double ddn = Dot(ray.dir, orienting_normal);
			const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

			if (cos2t < 0.0) { // �S���˂���
				return direct_light + Multiply(transmittance_ratio, Multiply(obj.mat.ref, (radiance(reflection_ray, medium, rng, depth + 1, maxDepth)))) / (1.0 - scattering_probability) / russian_roulette_probability;
			}
			// ���܂��Ă�������
			Vec tdir = Normalize(ray.dir * nnt - normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(cos2t)));

			// Schlick�ɂ��Fresnel�̔��ˌW���̋ߎ�
			const double a = nt - nc, b = nt + nc;
			const double R0 = (a * a) / (b * b);
			const double c = 1.0 - (into ? -ddn : Dot(tdir, normal));
			const double Re = R0 + (1.0 - R0) * pow(c, 5.0);
			const double Tr = 1.0 - Re; // ���܌��̉^�Ԍ��̗�
			const double probability = 0.25 + 0.5 * Re;

			// ���ܕ�������̒��ڌ��T���v�����O����
			Ray refraction_ray = Ray(hitpoint, tdir);
			intersect_scene_triangle(refraction_ray, &lintersect);
			Vec direct_light_refraction;
			if (lintersect.obj_id == LightID) {
				direct_light_refraction = direct_radiance(hitpoint, -1.0 * orienting_normal, obj_id, tri_num, refraction_ray.org + lintersect.hitpoint.distance * refraction_ray.dir, lintersect);
			}

			// ���ȏヌ�C��ǐՂ�������܂Ɣ��˂̂ǂ��炩�����ǐՂ���B�i�����Ȃ��Ǝw���I�Ƀ��C��������j
			// ���V�A�����[���b�g�Ō��肷��B
			if (rng.next01() < probability) { // ����
				return direct_light +
					Multiply(transmittance_ratio, Multiply(obj.mat.ref, radiance(reflection_ray, medium, rng, depth + 1, maxDepth) * Re))
					/ probability
					/ (1.0 - scattering_probability)
					/ russian_roulette_probability;
			}
			else { // ����
				//   // �������������̂Ȃ�medium���ω�
				//Medium next_medium = medium;
				//if (obj.mat.type == TRANSLUCENT && into) {
				//	next_medium = triangles[obj_id][tri_num].mat.medium;
				//}
				//else if (obj.mat.type == TRANSLUCENT && !into) {
				//	next_medium = Medium();
				//}

				return direct_light_refraction +
					Multiply(transmittance_ratio, Multiply(obj.mat.ref, radiance(Ray(hitpoint, tdir), medium, rng, depth + 1, maxDepth) * Tr))
					/ (1.0 - probability)
					/ (1.0 - scattering_probability)
					/ russian_roulette_probability;
			}}break;
		case TRANSLUCENT: {
			Ray reflection_ray = Ray(hitpoint, ray.dir - normal * 2.0 * Dot(normal, ray.dir));

			// ���˕�������̒��ڌ��T���v�����O����
			Intersection lintersect;
			intersect_scene_triangle(reflection_ray, &lintersect);
			Vec direct_light;
			if (lintersect.obj_id == LightID) {
				direct_light = direct_radiance(hitpoint, orienting_normal, obj_id, tri_num, reflection_ray.org + lintersect.hitpoint.distance * reflection_ray.dir, lintersect);
			}
			bool into = Dot(normal, orienting_normal) > 0.0; // ���C���I�u�W�F�N�g����o��̂��A����̂��A�o��Ȃ��false

															 // Snell�̖@��
			const double nc = 1.0; // �^��̋��ܗ�
			const double nt = 1.3; // �I�u�W�F�N�g�̋��ܗ�
			const double nnt = into ? nc / nt : nt / nc;
			const double ddn = Dot(ray.dir, orienting_normal);
			const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

			if (cos2t < 0.0) { // �S���˂���
				return direct_light + Multiply(transmittance_ratio, Multiply(obj.mat.ref, (radiance(reflection_ray, medium, rng, depth + 1, maxDepth)))) / (1.0 - scattering_probability) / russian_roulette_probability;
			}
			// ���܂��Ă�������
			Vec tdir = Normalize(ray.dir * nnt - normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(cos2t)));

			// Schlick�ɂ��Fresnel�̔��ˌW���̋ߎ�
			const double a = nt - nc, b = nt + nc;
			const double R0 = (a * a) / (b * b);
			const double c = 1.0 - (into ? -ddn : Dot(tdir, normal));
			const double Re = R0 + (1.0 - R0) * pow(c, 5.0);
			const double Tr = 1.0 - Re; // ���܌��̉^�Ԍ��̗�
			const double probability = 0.25 + 0.5 * Re;

			// ���ܕ�������̒��ڌ��T���v�����O����
			Ray refraction_ray = Ray(hitpoint, tdir);
			intersect_scene_triangle(refraction_ray, &lintersect);
			Vec direct_light_refraction;
			if (lintersect.obj_id== LightID) {
				direct_light_refraction = direct_radiance(hitpoint, -1.0 * orienting_normal, obj_id,tri_num, refraction_ray.org + lintersect.hitpoint.distance * refraction_ray.dir,lintersect);
			}

			// ���ȏヌ�C��ǐՂ�������܂Ɣ��˂̂ǂ��炩�����ǐՂ���B�i�����Ȃ��Ǝw���I�Ƀ��C��������j
			// ���V�A�����[���b�g�Ō��肷��B
			if (rng.next01() < probability) { // ����
				return direct_light +
					Multiply(transmittance_ratio, Multiply(obj.mat.ref, radiance(reflection_ray, medium, rng, depth + 1, maxDepth) * Re))
					/ probability
					/ (1.0 - scattering_probability)
					/ russian_roulette_probability;
			}
			else { // ����
				   // �������������̂Ȃ�medium���ω�
				Medium next_medium = medium;
				if (obj.mat.type == TRANSLUCENT && into) {
					next_medium = triangles[obj_id][tri_num].mat.medium;
				}
				else if (obj.mat.type == TRANSLUCENT && !into) {
					next_medium = Medium();
				}

				return direct_light_refraction +
					Multiply(transmittance_ratio, Multiply(obj.mat.ref, radiance(Ray(hitpoint, tdir), next_medium, rng, depth + 1, maxDepth) * Tr))
					/ (1.0 - probability)
					/ (1.0 - scattering_probability)
					/ russian_roulette_probability;
			}
		} break;
	}
	}

	Assertion(false, "Error!!");
}

inline double clamp(double x) {
	if (x < 0.0)
		return 0.0;
	if (x > 1.0)
		return 1.0;
	return x;
}

inline int to_int(double x) {
	return int(pow(clamp(x), 1 / 2.2) * 255 + 0.5);
}

// PPM�t�@�C���̕ۑ�
void save_ppm_file(const std::string &filename, const Color *image, const int width, const int height) {
	std::ofstream writer(filename.c_str(), std::ios::out);
	writer << "P3" << std::endl;
	writer << width << " " << height << std::endl;
	writer << 255 << std::endl;
	for (int i = 0; i < width * height; i++) {
		const int r = to_int(image[i].x);
		const int g = to_int(image[i].y);
		const int b = to_int(image[i].z);
		writer << r << " " << g << " " << b << " ";
	}
	writer.close();
}

void save_box_obj_blender_file(const std::string &filename,Vec *v) {
	std::ofstream writer(filename.c_str(), std::ios::out);
	writer << "# Blender v2.78 (sub 0) OBJ File: ''# www.blender.orgo Cube\n";

	for (int i = 0; i <8 ; i++) {
		writer << "v" << " " << std::fixed << std::setprecision(5) << (float)(v[i].x) << " " << std::fixed << std::setprecision(5) << (float)(v[i].y) << " " << std::fixed << std::setprecision(5) << (float)(v[i].z) << std::endl;
	}
	
	writer << "vn -1.0000 0.0000 0.0000" << std::endl;
	writer << "vn 0.0000 0.0000 -1.0000" << std::endl;
		writer<< "vn 1.0000 0.0000 0.0000" << std::endl; 
		writer<<"vn 0.0000 0.0000 1.0000" << std::endl; 
		writer << "vn 0.0000 -1.0000 0.0000"<< std::endl; 
		writer << "vn 0.0000 1.0000 0.0000"<<std::endl;
		writer << "s off" << std::endl;
		writer << "f 2//1 3//1 1//1" << std::endl;
		writer << "f 4//2 7//2 3//2" << std::endl;
		writer << "f 8//3 5//3 7//3" << std::endl;
		writer << "f 6//4 1//4 5//4" << std::endl;
		writer << "f 7//5 1//5 3//5" << std::endl;
		writer << "f 4//6 6//6 8//6" << std::endl;
		writer << "f 2//1 4//1 3//1" << std::endl;
		writer << "f 4//2 8//2 7//2" << std::endl;
		writer << "f 8//3 6//3 5//3" << std::endl;
		writer << "f 6//4 2//4 1//4" << std::endl;
		writer << "f 7//5 5//5 1//5" << std::endl;
		writer << "f 4//6 2//6 6//6" ;
	writer.close();
}

// �i�s�x��\�����郁�\�b�h
inline void progressBar(int x, int total, int width = 50) {
	double ratio = (double)x / total;
	int tick = (int)(width * ratio);
	std::string bar(width, ' ');
	std::fill(bar.begin(), bar.begin() + tick, '+');
	printf("[ %6.2f %% ] [ %s ]", 100.0 * ratio, bar.c_str());
	printf("%c", x >= total ? '\n' : '\r');
}



//object��load����֐�


//�K��light��0�Ֆڂɓǂݍ���
inline void objectload(int i , std::vector<std::string> strList) {
	
	std::string inputfile = strList[i];
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;

	std::string err;
	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, inputfile.c_str());

	if (!err.empty()) { // `err` may contain warning message.
		std::cerr << err << std::endl;
	}

	if (!ret) {
		exit(1);
	}

	// Loop over shapes
	for (size_t s = 0; s < shapes.size(); s++) {
		// Loop over faces(polygon)
		int  face_number = shapes[s].mesh.num_face_vertices.size();
		triangles[i].resize(face_number);

		size_t index_offset = 0;

		for (size_t f = 0; f < face_number; f++) {
			int fv = shapes[s].mesh.num_face_vertices[f];
			// Loop over vertices in the face.
			for (size_t ve = 0; ve < fv; ve++) {
				// access to vertex

				tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + ve];
				float vx = attrib.vertices[3 * idx.vertex_index + 0];
				float vy = attrib.vertices[3 * idx.vertex_index + 1];
				float vz = attrib.vertices[3 * idx.vertex_index + 2];
				float nx = attrib.normals[3 * idx.normal_index + 0];
				float ny = attrib.normals[3 * idx.normal_index + 1];
				float nz = attrib.normals[3 * idx.normal_index + 2];
				/*float tx = attrib.texcoords[2 * idx.texcoord_index + 0];
				float ty = attrib.texcoords[2 * idx.texcoord_index + 1];*///���ꂪ����ƂȂ񂩃G���[���o��

				Vec  vertex = Vec(vx, vy, vz);
				Vec  normal = Vec(nx, ny, nz);

				triangles[i][f].v[ve] = vertex;
				triangles[i][f].n[ve] = normal;
				

							if (i == 0) {
								triangles[i][f].mat = blueMat;//TODO �e�I�u�W�F�N�g�ɂ��ĕύX�ł���悤�ɂ���
							}
							//else if(i==1){
							//	triangles[i][f].mat =redMat;//TODO �e�I�u�W�F�N�g�ɂ��ĕύX�ł���悤�ɂ���
							//}
							//else if (i == 2) {
							//	triangles[i][f].mat = blueMat;//TODO �e�I�u�W�F�N�g�ɂ��ĕύX�ł���悤�ɂ���
							//}
							//else if (i == 3) {
							//	triangles[i][f].mat = grayMat;//TODO �e�I�u�W�F�N�g�ɂ��ĕύX�ł���悤�ɂ���
							//}
							//else if (i == 4) {
							//	triangles[i][f].mat = grayMat;//TODO �e�I�u�W�F�N�g�ɂ��ĕύX�ł���悤�ɂ���
							//}
							//else if (i == 5) {
							//	triangles[i][f].mat = grayMat;//TODO �e�I�u�W�F�N�g�ɂ��ĕύX�ł���悤�ɂ���
							//}
							//else if (i == 6) {
							//	triangles[i][f].mat = diamondMat;//TODO �e�I�u�W�F�N�g�ɂ��ĕύX�ł���悤�ɂ���
							//}
			}

			triangles[i][f].bbox[0][0] = std::min(std::min(triangles[i][f].v[0].x, triangles[i][f].v[1].x), triangles[i][f].v[2].x);
			triangles[i][f].bbox[0][1] = std::min(std::min(triangles[i][f].v[0].y, triangles[i][f].v[1].y), triangles[i][f].v[2].y);
			triangles[i][f].bbox[0][2] = std::min(std::min(triangles[i][f].v[0].z, triangles[i][f].v[1].z), triangles[i][f].v[2].z);
			triangles[i][f].bbox[1][0] = std::max(std::max(triangles[i][f].v[0].x, triangles[i][f].v[1].x), triangles[i][f].v[2].x);
			triangles[i][f].bbox[1][1] = std::max(std::max(triangles[i][f].v[0].y, triangles[i][f].v[1].y), triangles[i][f].v[2].y);
			triangles[i][f].bbox[1][2] = std::max(std::max(triangles[i][f].v[0].z, triangles[i][f].v[1].z), triangles[i][f].v[2].z);
			triangles[i][f].normal = Normalize(Cross((triangles[i][f].v[1]- triangles[i][f].v[0]), (triangles[i][f].v[2] - triangles[i][f].v[0])));
			index_offset += fv;

			// per-face material
			shapes[s].mesh.material_ids[f];
		}
	}

	std::cout<< objList[i] << "���ǂݍ��܂�܂����B"<<std::endl;

}









// ���C���֐�
int main(int argc, char **argv) {
	
	triangles.resize(objList.size());

	//scene
	for (int i = 0; i < objList.size(); i++) {
		objectload(i, objList);
	}
	
	

	// �R�}���h�����̃p�[�X
	int width = 640;
	int height = 480;
	int samples =  1;
	int maxDepth = 1;
	/*for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "--width") == 0) {
			width = std::atoi(argv[++i]);
		}

		if (strcmp(argv[i], "--height") == 0) {
			height = std::atoi(argv[++i]);
		}

		if (strcmp(argv[i], "--samples") == 0) {
			samples = std::atoi(argv[++i]);
		}

		if (strcmp(argv[i], "--depth") == 0) {
			maxDepth = std::atoi(argv[++i]);
		}
	}
*/
	// �p�����[�^�̕\��
	printf("-- Parameters --\n");
	printf("    Width: %d\n", width);
	printf("   Height: %d\n", height);
	printf("  Samples: %d\n", samples);
	printf("Max depth: %d\n", maxDepth);
	printf("\n");

	// �J����
	const Vec camPos = 0.1 * Vec(49.0, 60.0, 295.6);
	const Vec camDir = Normalize(Vec(-0.045, -0.042612, -1.0));
	const Ray camera(camPos, camDir);

	// �X�N���[���̊��x�N�g��
	const Vec cx = Vec(width * 0.5135 / height, 0.0, 0.0);
	const Vec cy = Normalize(Cross(cx, camera.dir)) * 0.5135;

	// �����_�����O���[�v
	auto image = std::make_unique<Color[]>(width * height);
	std::atomic<int> progress(0);
	std::mutex mtx;
	//parallel_for(0, height, [&](int y) {
	//	for (int x = 0; x < width; x++) {
	//		// ����
	//		Random rng(y * width + x);

	//		// �s�N�Z���F�̏�����
	//		Color &pixel = image[y * width + x];
	//		pixel = Color(0.0, 0.0, 0.0);

	//		// �T���v�����O
	//		for (int s = 0; s < samples; s++) {
	//			// �e���g�t�B���^�[�ɂ���ăT���v�����O
	//			const double r1 = 2.0 * rng.next01();
	//			const double r2 = 2.0 * rng.next01();
	//			const double dx = r1 < 1.0 ? sqrt(r1) - 1.0 : 1.0 - sqrt(2.0 - r1);
	//			const double dy = r2 < 1.0 ? sqrt(r2) - 1.0 : 1.0 - sqrt(2.0 - r2);

	//			const double px = (x + dx + 0.5) / width - 0.5;
	//			const double py = ((height - y - 1) + dy + 0.5) / height - 0.5;

	//			// ���ˋP�x�̌v�Z
	//			const Vec dir = cx * px + cy * py + camera.dir;
	//			const Ray ray(camera.org + dir * 13.0, Normalize(dir));
	//			const Color L = radiance(ray, Medium(), rng, 0, maxDepth);
	//			Assertion(L.isValid(), "Radiance is invalid: (%f, %f %f)", L.x, L.y, L.z);

	//			pixel = pixel + L;
	//		}
	//		pixel = pixel / samples;

	//		// �i�s�x�̕\��
	//		mtx.lock();
	//		progressBar(++progress, width * height);
	//		mtx.unlock();
	//	}
	//});

	// PPM�t�@�C����ۑ�
	//const std::string outfile = std::string(OUTPUT_DIRECTORY) + "image.ppm";
	//save_ppm_file(outfile, image.get(), width, height);

	 std::vector<TRIANGLE*> polygons;
	for (int i = 0; i < objList.size(); i++) {
		for (int j = 0; j < triangles[i].size(); j++) {
			polygons.push_back(&triangles[i][j]);
		}
	}

	Vec v[8];

	constructBVH(polygons);

	//nodes�̊e�X��bbox��obj�t�@�C���œf���o��(blender�ł����ɂ�)
	for (int i = 0; i < (sizeof(nodes) / (sizeof(BVH_node)));i++) {
		
		/*v[2] = Vec(1.617, 1.2721, 7.44);
		v[5] = Vec(6.413,6.022,11.165 );*/

		v[2] = Vec(nodes[i].bbox[0][0], nodes[i].bbox[0][1], nodes[i].bbox[0][2]);
		v[5] = Vec(nodes[i].bbox[1][0], nodes[i].bbox[1][1], nodes[i].bbox[1][2]);
		v[0] = v[2] + Vec (0.0, 0.0, v[5].z-v[2].z);
		v[1] = v[2] + Vec(0.0, v[5].y - v[2].y, v[5].z - v[2].z);
		v[3] = v[2] + Vec(0.0, v[5].y - v[2].y, 0.0);
		v[4] = v[2] + Vec(v[5].x-v[2].x, 0.0, v[5].z - v[2].z);
		v[6] = v[2] + Vec(v[5].x - v[2].x, 0.0, 0.0);
		v[7] = v[2] + Vec(v[5].x - v[2].x, v[5].y - v[2].y, 0.0);
		std::string fname = "bbox";
		std::ostringstream ss;
		ss << i;
		fname = fname + ss.str();
		fname = fname + ".obj";
		const std::string outfile = std::string(OUTPUT_DIRECTORY) + fname;
		save_box_obj_blender_file(outfile, v);
	

	}
}
