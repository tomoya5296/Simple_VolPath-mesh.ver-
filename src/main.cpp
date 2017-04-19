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
#include<time.h>
#include <fenv.h>
#include "common.h"
#include "dirs.h"
#include "parallel.h"
#include "vec.h"
#include "ray.h"
#include "sphere.h"
#include "random.h"
#include "TRIANGLE.h"
#include"scene.h"
#include"tiny_obj_loader.h"
#include "bvh.h"






// シーンの補助情報
const int LightID = 0;
std::vector<std::string> objList =
{
	"light_plane.obj","left_plane.obj","right_plane.obj","up_plane.obj","bottom_plane.obj","far_plane.obj","bunny.obj"
};
std::vector<std::vector <TRIANGLE>>triangles;//三角形たちを二次元配列で持つ


int Intersect(const std::vector<TRIANGLE *> &polygons, const Ray ray, Intersection  *intersection) {
	for (int i = 0; i<polygons.size(); i++) {
		TRIANGLE *polygon = polygons[i];
		Hitpoint hitpoint;
		if (polygon->intersect(ray, &hitpoint)) {
			if (hitpoint.distance <= intersection->hitpoint.distance) {
				intersection->hitpoint = hitpoint;
				intersection->Mat = hitpoint.tri->mat;
				intersection->hitpoint.orienting_normal = Dot(intersection->hitpoint.normal, ray.dir) < 0.0 ? intersection->hitpoint.normal : (-1.0 * intersection->hitpoint.normal);
			}
		}
	}
	return -1;
}



bool IntersectAABB(float aabb[2][3], const Ray &ray) {
	double ray_dir[3], ray_org[3];

	for (int i = 0; i < 3; i++) {
		if (i == 0) {
			ray_dir[i] = ray.dir.x;
			ray_org[i] = ray.org.x;
		}
		else if (i == 1) {
			ray_dir[i] = ray.dir.y;
			ray_org[i] = ray.org.y;
		}
		else {
			ray_dir[i] = ray.dir.z;
			ray_org[i] = ray.org.z;
		}
	}

	double t_max = INF;
	double t_min = -INF;

	for (int i = 0; i < 3; i++) 
	{
		
		double t1 = (aabb[0][i] - ray_org[i]) /( ray_dir[i]+EPS);
feclearexcept(FE_ALL_EXCEPT);

		double t2 = (aabb[1][i] - ray_org[i]) /( ray_dir[i]+EPS);
		double t_far = std::max(t1, t2);
		double t_near = std::min(t1, t2);
		t_max = std::min(t_max, t_far);
		t_min = std::max(t_min, t_near);
		if (t_min > t_max) { return false; }
	}
	return true;
};




TRIANGLE * Intersect(BVH_node *nodes, int index, const Ray &ray, Intersection *intersect) {
	// AABB とレイの交差判定
	if (IntersectAABB(nodes[index].bbox, ray)) {
		// 交差している

		// 中間ノードか？
		if (nodes[index].children[0] != -1) {
			// 中間ノード
			// 両方の子ノードについて調べる
			TRIANGLE *childResult = nullptr;
			for (int i = 0; i<2; i++) {
				TRIANGLE *result = Intersect(nodes, nodes[index].children[i],  ray, intersect);
				if (result != nullptr) {
					childResult = result;
				}
			}
			if (childResult != nullptr) return childResult;
		}
		else {
			// 葉ノード
			TRIANGLE *result = nullptr;
			for (TRIANGLE *tri : nodes[index].polygons) {
				// ポリゴンとレイの交差判定
				// distance に交差していた場合のレイからの距離、normal にポリゴンの法線が入る
				Hitpoint hitopoint;
				if (tri->intersect(ray,&hitopoint)) {
					// 既に交差したと判定された他のポリゴンより、レイの始点に近いかどうか
					if ( hitopoint.distance<intersect->hitpoint.distance ) {
						result = tri;
						intersect->hitpoint=hitopoint;
						intersect->Mat = tri->mat;
						intersect->obj_id = tri->obj_id;
					}
				}
			}
			if (result != nullptr) return result;
		}
	}
	else {
		return nullptr;// 交差していない (何もする必要なし)
	}
	return nullptr;
}



//光の向きが一意(SPECULAR,REFRACTION)のとき使う//注意　光源のどこかをinteersectで必ず渡すこと:取れていなければ使わない ref)SPECULAR

Color direct_radiance(const Vec &v0, const Vec &normal, const TRIANGLE* tri, const Vec &light_pos,Intersection lintersect) {
	//tri は光源ではなくもとのobjのメッシュ
	const Vec light_normal = lintersect.hitpoint.tri->normal;
	const Vec light_dir = Normalize(light_pos - v0);
	const double dist2 = (light_pos - v0).LengthSquared();
	const double dot0 = Dot(normal, light_dir);
	const double dot1 = Dot(light_normal, -1.0 * light_dir);

	if (dot0 >= 0.0 && dot1 >= 0.0) {
		const double G = dot0 * dot1 / dist2;
		TRIANGLE* HitTri = nullptr;
		Intersection intersect;
		HitTri = Intersect(nodes, 0, Ray(v0, light_dir), &intersect);
		if (std::abs(sqrt(dist2) - intersect.hitpoint.distance) < 1e-3) {
			Vec edge1 = (intersect.hitpoint.tri->v[1]- intersect.hitpoint.tri->v[0]);
			Vec edge2 = (intersect.hitpoint.tri->v[2] - intersect.hitpoint.tri->v[0]);
			double S = Cross(edge1,edge2).Length()/2;//lightのメッシュ一つの表面積
			return Multiply(tri->mat.ref, intersect.hitpoint.tri->mat.Le) * (1.0 / PI) * G / (1.0 / (triangles[LightID].size()*S));
			//(1.0 / PI)てなんでだっけ
		}
	}
	return Color(0.0);
}



// 光源上の点をサンプリングして直接光を計算する。//DIFFUSE面で用いる
Color direct_radiance_sample(const Vec &v0, const Vec &normal, const TRIANGLE* tri, double u0, double u1,double u2) {
	// 光源上の一点をサンプリングする
	//u0,u1は0～1の乱数
	int r1 = (int)(u0*(2));//TODO メッシュの数を撃ち込まなくてもできるようにする.
	TRIANGLE *light_triangle =&triangles[LightID][r1];//ランダムにlightのメッシュが取れた
	Vec light_pos = light_triangle->v[0] + u1*(light_triangle->v[1] - light_triangle->v[0]) + u2*(1.0-u1)*(light_triangle->v[2] - light_triangle->v[1]);
	Vec dir = Normalize(light_pos - v0);

	TRIANGLE* HitTri = nullptr;
	Intersection lintersect;
	HitTri = Intersect(nodes, 0, Ray(v0, dir), &lintersect);


	if (HitTri != nullptr) {
		if (HitTri->obj_id == LightID) {
			return direct_radiance(v0, normal, tri, light_pos, lintersect);
		}
	}
	return Color(0.0);
}

//// 指定した位置への直接光を計算する。ただし、空間中の点が光を受ける。
Color direct_radiance_media(const Vec &v0, const Vec &light_pos,const TRIANGLE &light_triangle) {
	const Vec light_normal = light_triangle.normal;
	const Vec light_dir = Normalize(light_pos - v0);
	const double dist2 = (light_pos - v0).LengthSquared();
	const double dot1 = Dot(light_normal, -1.0 * light_dir);

	if (dot1 >= 0) {

		Vec edge1 = (light_triangle.v[1] - light_triangle.v[0]);
		Vec edge2 = (light_triangle.v[2] - light_triangle.v[0]);
		double S = Cross(edge1, edge2).Length() / 2;//lightのメッシュ一つの表面積
		const double G = dot1 / dist2;
		TRIANGLE* HitTri = nullptr;
		Intersection lintersect;
		HitTri = Intersect(nodes, 0, Ray(v0, light_dir), &lintersect);
		if (fabs(sqrt(dist2) - lintersect.hitpoint.distance) < 1e-3) {
			const Vec ret = light_triangle.mat.Le * (1.0 / PI) * G / (1.0 / 2*S);
			return ret;
		}
	}
	return Color(0.0);
}

// 媒質内にある点v0に対する直接光の影響を計算する
Color direct_radiance_sample_media(const Vec &v0, const Medium &mat, double u0, double u1,double u2) {
	// 光源上の一点をサンプリングする
	int r1 = (int)(u0*(2));//TODO メッシュの数を撃ち込まなくてもできるようにする.
	TRIANGLE light_triangle =triangles[LightID][r1];//ランダムにlightのメッシュが取れた
	Vec light_pos = light_triangle.v[0] + u1*(light_triangle.v[1] - light_triangle.v[0]) + u2*(1.0-u1)*(light_triangle.v[2] - light_triangle.v[1]);
	Vec dir = Normalize(light_pos - v0);
	TRIANGLE* HitTri = nullptr;
	Intersection lintersect;
	HitTri = Intersect(nodes, 0, Ray(v0,dir), &lintersect);
	// 媒質の外に出るまでの光の減衰を計算
	const Color sigT = mat.sigS + mat.sigA;
	const Vec transmittance_ratio = Vec::exp(-sigT * lintersect.hitpoint.distance);
	// 外に出る点への直接光の影響を計算する
	const Vec v1 = v0 + lintersect.hitpoint.distance * dir;
	const Color direct_light = direct_radiance_media(v1, light_pos,light_triangle);
	return Multiply(transmittance_ratio, direct_light);
}



// ray方向からの放射輝度を求める
// ボリュームレンダリング方程式に基づく
Color radiance(const Ray &ray, const Medium &medium, Random &rng, int depth, int maxDepth) {
	TRIANGLE* HitTri = nullptr;
	Intersection intersection;
	HitTri = Intersect(nodes, 0, ray, &intersection);
	if (HitTri == nullptr) {
		return Color(0.0);
	}
	const double t = intersection.hitpoint.distance;
	const TRIANGLE &obj = *HitTri;
	const Vec normal = obj.normal;//ok
	const Vec orienting_normal = Dot(normal, ray.dir) < 0.0 ? normal : (-1.0 * normal); // 交差位置の法線（物体からのレイの入出を考慮）
	const Vec hitpoint = ray.org +t*ray.dir;//ok
	// 最大のbounce数を評価
	if (depth >= maxDepth) {
		if (Dot(-ray.dir,normal)>0.0) {
			return obj.mat.Le;
		}
		else {
			return Color(0.0);
		}
	}

	
	// 一定以上レイを追跡したらロシアンルーレットを実行し追跡を打ち切るかどうかを判断する
	double russian_roulette_probability = std::max(obj.mat.ref.x, std::max(obj.mat.ref.y, obj.mat.ref.z));
	russian_roulette_probability = std::max(0.05, russian_roulette_probability);
	if (depth > 5) {
		if (rng.next01() > russian_roulette_probability) {
			if (Dot(-ray.dir, normal) > 0.0) {
				return obj.mat.Le / (1.0 - russian_roulette_probability);
			}
			else {
				return Color(0.0);
			}
		}
	}
	else {
		russian_roulette_probability = 1.0;
	}
	// レイと物体の交差点からの放射輝度を計算するか、途中の点における周囲からの影響を計算するかを選択するロシアンルーレット
	const Color sigT = medium.sigS + medium.sigA;
	const double tr_average = (sigT.x + sigT.y + sigT.z) / 3.0;
	const double sc_average = (medium.sigS.x + medium.sigS.y + medium.sigS.z) / 3.0;
	double scattering_probability = std::max(0.0, std::min(sc_average, 0.90));
	if (sigT.isZero()) {
		scattering_probability = 0.0;
	}

	if (rng.next01() < scattering_probability) {
		//平均自由工程
		const double u = rng.next01();
		const double d = -log(1.0 + u * (exp(-tr_average * t) - 1.0)) / tr_average;
		const double pdf = exp(-tr_average * d) * (-tr_average / (exp(-tr_average * t) - 1.0));

		//次レイの方向を球面上から一様サンプリング
		double r1 = 2.0 * PI * rng.next01();
		double r2 = 1.0 - 2.0 * rng.next01();
		const Vec next_dir(sqrt(1.0 - r2*r2) * cos(r1), sqrt(1.0 - r2*r2) * sin(r1), r2);
		const Ray next_ray(ray.org + d * ray.dir, next_dir);

		const Vec transmittance_ratio = Vec::exp(-sigT * d);
		const Vec direct_light = direct_radiance_sample_media(next_ray.org, medium, rng.next01(), rng.next01(),rng.next01());

		// 位相関数
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
		// レイと物体の交差点からの放射輝度伝達を計算
		const Vec transmittance_ratio = Vec::exp(-sigT * t);
		switch (obj.mat.type) {
		case DIFFUSE: {
			// 直接光のサンプリングを行う
			if (obj.obj_id!=LightID) {
				const int shadow_ray = 1;
				Vec direct_light;
				for (int i = 0; i < shadow_ray; i++) {
					direct_light = direct_light + direct_radiance_sample(hitpoint, orienting_normal, &obj, rng.next01(), rng.next01(),rng.next01()) / shadow_ray;
				}
				
				// orienting_normalの方向を基準とした正規直交基底(w, u, v)を作る。この基底に対する半球内で次のレイを飛ばす。
				Vec w, u, v;
				w = orienting_normal;
				if (std::abs(w.x) > 0.1)
					u = Normalize(Cross(Vec(0.0, 1.0, 0.0), w));
				else
					u = Normalize(Cross(Vec(1.0, 0.0, 0.0), w));
				v = Cross(w, u);

				// コサイン項を使った重点的サンプリング
				const double r1 = 2.0 * PI * rng.next01();
				const double r2 = rng.next01();
				const double r2s = sqrt(r2);
				Vec dir = Normalize((u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1.0 - r2)));
				// 減衰を考慮しつつ次を計算
				return direct_light + Multiply(transmittance_ratio, Multiply(obj.mat.ref, radiance(Ray(hitpoint, dir), medium, rng, depth + 1, maxDepth))) / (1.0 - scattering_probability) / russian_roulette_probability;
			}
			else if (depth == 0) {
				if (Dot(-ray.dir,obj.normal)<0.0) {
					return obj.mat.Le;
				}
				else {
					return Color(0.0);
				}
			}
			else {
				return Color(0.0);
			}
		} break;

		case SPECULAR: {
			// 完全鏡面なのでレイの反射方向は決定的。
			// ロシアンルーレットの確率で除算するのは上と同じ。
			//Intersection lintersect;//反射光の情報
			Ray reflection_ray = Ray(hitpoint+Vec(EPS), Normalize(ray.dir + normal * 2.0 * Dot(-normal, ray.dir)));
			return	Multiply(transmittance_ratio, radiance(reflection_ray, medium, rng, depth + 1, maxDepth)) / (1.0 - scattering_probability) / russian_roulette_probability;
		} break;
		case REFRACTION: {			
			Ray reflection_ray = Ray(hitpoint, ray.dir - normal * 2.0 * Dot(normal, ray.dir));
			// 反射方向からの直接光サンプリングする
			Intersection llintersect;
			TRIANGLE* HitTri = nullptr;
			HitTri = Intersect(nodes, 0, reflection_ray,
				&llintersect);
			//intersect_scene_triangle(reflection_ray, &lintersect);
			Vec direct_light;
			if (llintersect.obj_id == LightID) {
				direct_light = direct_radiance(hitpoint, orienting_normal,HitTri, reflection_ray.org + llintersect.hitpoint.distance * reflection_ray.dir, llintersect);
			}
			bool into = Dot(normal, orienting_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか、出るならばfalse

															 // Snellの法則
			const double nc = 1.0; // 真空の屈折率
			const double nt = 1.3; // オブジェクトの屈折率
			const double nnt = into ? nc / nt : nt / nc;
			const double ddn = Dot(ray.dir, orienting_normal);
			const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
			if (cos2t < 0.0) { // 全反射した
				return direct_light + Multiply(transmittance_ratio, Multiply(obj.mat.ref, (radiance(reflection_ray, medium, rng, depth + 1, maxDepth)))) / (1.0 - scattering_probability) / russian_roulette_probability;
			}
			// 屈折していく方向
			Vec tdir = Normalize(ray.dir * nnt - normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(cos2t)));
			// SchlickによるFresnelの反射係数の近似
			const double a = nt - nc, b = nt + nc;
			const double R0 = (a * a) / (b * b);
			const double c = 1.0 - (into ? -ddn : Dot(tdir, normal));
			const double Re = R0 + (1.0 - R0) * pow(c, 5.0);
			const double Tr = 1.0 - Re; // 屈折光の運ぶ光の量
			const double probability = 0.25 + 0.5 * Re;

			// 屈折方向からの直接光サンプリングする
			Ray refraction_ray = Ray(hitpoint, tdir);


			Intersection lintersect;
			HitTri = nullptr;
			HitTri = Intersect(nodes, 0, reflection_ray, &lintersect);
			//intersect_scene_triangle(refraction_ray, &lintersect);
			Vec direct_light_refraction;


			if (lintersect.obj_id == LightID) {
				direct_light_refraction = direct_radiance(hitpoint, -1.0 * orienting_normal, HitTri, refraction_ray.org + lintersect.hitpoint.distance * refraction_ray.dir, lintersect);
			}

			// 一定以上レイを追跡したら屈折と反射のどちらか一方を追跡する。（さもないと指数的にレイが増える）
			// ロシアンルーレットで決定する。
			if (rng.next01() < probability) { // 反射
				return direct_light +
					Multiply(transmittance_ratio, Multiply(obj.mat.ref, radiance(reflection_ray, medium, rng, depth + 1, maxDepth) * Re))
					/ probability
					/ (1.0 - scattering_probability)
					/ russian_roulette_probability;
			}
			else { // 屈折
				return direct_light_refraction +
					Multiply(transmittance_ratio, Multiply(obj.mat.ref, radiance(Ray(hitpoint, tdir), medium, rng, depth + 1, maxDepth) * Tr))
					/ (1.0 - probability)
					/ (1.0 - scattering_probability)
					/ russian_roulette_probability;
			}}break;
		case TRANSLUCENT: {
			Ray reflection_ray = Ray(hitpoint, ray.dir - normal * 2.0 * Dot(normal, ray.dir));

			// 反射方向からの直接光サンプリングする
			Intersection lintersect;
			TRIANGLE* HitTri = nullptr;
			HitTri = Intersect(nodes, 0, reflection_ray, &lintersect);
			Vec direct_light;
			if (lintersect.obj_id == LightID) {
				direct_light = direct_radiance(hitpoint, orienting_normal,HitTri, reflection_ray.org + lintersect.hitpoint.distance * reflection_ray.dir, lintersect);
			}
			bool into = Dot(normal, orienting_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか、出るならばfalse

															 // Snellの法則
			const double nc = 1.0; // 真空の屈折率
			const double nt = 1.3; // オブジェクトの屈折率
			const double nnt = into ? nc / nt : nt / nc;
			const double ddn = Dot(ray.dir, orienting_normal);
			const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

			if (cos2t < 0.0) { // 全反射した
				return direct_light + Multiply(transmittance_ratio, Multiply(obj.mat.ref, (radiance(reflection_ray, medium, rng, depth + 1, maxDepth)))) / (1.0 - scattering_probability) / russian_roulette_probability;
			}
			// 屈折していく方向
			Vec tdir = Normalize(ray.dir * nnt - normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(cos2t)));

			// SchlickによるFresnelの反射係数の近似
			const double a = nt - nc, b = nt + nc;
			const double R0 = (a * a) / (b * b);
			const double c = 1.0 - (into ? -ddn : Dot(tdir, normal));
			const double Re = R0 + (1.0 - R0) * pow(c, 5.0);
			const double Tr = 1.0 - Re; // 屈折光の運ぶ光の量
			const double probability = 0.25 + 0.5 * Re;

			// 屈折方向からの直接光サンプリングする
			Ray refraction_ray = Ray(hitpoint, tdir);
			HitTri = Intersect(nodes, 0, reflection_ray, &lintersect);
			Vec direct_light_refraction;
			if (HitTri != nullptr) {
				if (lintersect.obj_id == LightID) {
					direct_light_refraction = direct_radiance(hitpoint, -1.0 * orienting_normal, HitTri, refraction_ray.org + lintersect.hitpoint.distance * refraction_ray.dir, lintersect);
				}
			}
			// 一定以上レイを追跡したら屈折と反射のどちらか一方を追跡する。（さもないと指数的にレイが増える）
			// ロシアンルーレットで決定する。
			if (rng.next01() < probability) { // 反射
				return direct_light +
					Multiply(transmittance_ratio, Multiply(obj.mat.ref, radiance(reflection_ray, medium, rng, depth + 1, maxDepth) * Re))
					/ probability
					/ (1.0 - scattering_probability)
					/ russian_roulette_probability;
			}
			else { // 屈折
				   // もし半透明物体ならmediumが変化
				Medium next_medium = medium;
				if (obj.mat.type == TRANSLUCENT && into) {
					next_medium = obj.mat.medium;
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

// PPMファイルの保存
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

void save_box_obj(const std::string &filename,Vec *v) {
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

// 進行度を表示するメソッド
inline void progressBar(int x, int total, int width = 50) {
	double ratio = (double)x / total;
	int tick = (int)(width * ratio);
	std::string bar(width, ' ');
	std::fill(bar.begin(), bar.begin() + tick, '+');
	printf("[ %6.2f %% ] [ %s ]", 100.0 * ratio, bar.c_str());
	printf("%c", x >= total ? '\n' : '\r');
}



//objectをloadする関数


//必ずlightは0盤目に読み込む
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

				Vec  vertex;
				Vec  normal;
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
				float ty = attrib.texcoords[2 * idx.texcoord_index + 1];*///これがあるとなんかエラーが出る
				vertex = Vec(vx, vy, vz);
				normal = Vec(nx, ny, nz);

				triangles[i][f].v[ve] = vertex;
				triangles[i][f].n[ve] = normal;
							if (i == 0) {
								triangles[i][f].mat = lightMat;//TODO 各オブジェクトについて変更できるようにする
							}
							else if(i==1){
								triangles[i][f].mat =redMat;//TODO 各オブジェクトについて変更できるようにする
							}
							else if (i == 2) {
								triangles[i][f].mat = blueMat;//TODO 各オブジェクトについて変更できるようにする
							}
							else if (i == 3) {
								triangles[i][f].mat = grayMat;//TODO 各オブジェクトについて変更できるようにする
							}
							else if (i == 4) {
								triangles[i][f].mat = grayMat;//TODO 各オブジェクトについて変更できるようにする
							}
							else if (i == 5) {
								triangles[i][f].mat = grayMat;//TODO 各オブジェクトについて変更できるようにする
							}
							else if (i == 6) {
								triangles[i][f].mat =blueMat;//TODO 各オブジェクトについて変更できるようにする
							}
			}

			triangles[i][f].bbox[0][0] = std::min(std::min(triangles[i][f].v[0].x, triangles[i][f].v[1].x), triangles[i][f].v[2].x);
			triangles[i][f].bbox[0][1] = std::min(std::min(triangles[i][f].v[0].y, triangles[i][f].v[1].y), triangles[i][f].v[2].y);
			triangles[i][f].bbox[0][2] = std::min(std::min(triangles[i][f].v[0].z, triangles[i][f].v[1].z), triangles[i][f].v[2].z);
			triangles[i][f].bbox[1][0] = std::max(std::max(triangles[i][f].v[0].x, triangles[i][f].v[1].x), triangles[i][f].v[2].x);
			triangles[i][f].bbox[1][1] = std::max(std::max(triangles[i][f].v[0].y, triangles[i][f].v[1].y), triangles[i][f].v[2].y);
			triangles[i][f].bbox[1][2] = std::max(std::max(triangles[i][f].v[0].z, triangles[i][f].v[1].z), triangles[i][f].v[2].z);
			triangles[i][f].normal =(triangles[i][f].n[0] + triangles[i][f].n[1] + triangles[i][f].n[2]) / 3; /*Normalize(Cross((triangles[i][f].v[1]- triangles[i][f].v[0]), (triangles[i][f].v[2] - triangles[i][f].v[0])));*/ 
			triangles[i][f].obj_id = i;
			index_offset += fv;

			// per-face material
			shapes[s].mesh.material_ids[f];
		}
	}

	std::cout<< objList[i] << "が読み込まれました。"<<std::endl;

}


// メイン関数
int main(int argc, char **argv) {
	
	triangles.resize(objList.size());

	//scene
	for (int i = 0; i < objList.size(); i++) {

		objectload(i, objList);

	}
	std::vector<TRIANGLE*> polygons;//BVHの関係で一次元配列で持たないといけない
	for (int i = 0; i < triangles.size(); i++) {
		for (int j = 0; j < triangles[i].size(); j++) {
			polygons.push_back(&triangles[i][j]);
		}
	}

	constructBVH(polygons);
	

	// コマンド引数のパース
	int width = 640;
	int height = 480;
	int samples = 10;
	int maxDepth = 50;
	for (int i = 1; i < argc; i++) {
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

	// パラメータの表示
	printf("-- Parameters --\n");
	printf("    Width: %d\n", width);
	printf("   Height: %d\n", height);
	printf("  Samples: %d\n", samples);
	printf("Max depth: %d\n", maxDepth);
	printf("\n");

	// カメラ
	const Vec camPos = 0.1 * Vec(49.0, 60.0, 295.6);
	const Vec camDir = Normalize(Vec(-0.045,-0.042612, -1.0));
	const Ray camera(camPos, camDir);

	// スクリーンの基底ベクトル
	const Vec cx = Vec(width * 0.5135 / height, 0.0, 0.0);
	const Vec cy = Normalize(Cross(cx, camera.dir)) * 0.5135;

	// レンダリングループ
	auto image = std::make_unique<Color[]>(width * height);
	std::atomic<int> progress(0);
	std::mutex mtx;

	clock_t start = clock();

	parallel_for(0, height, [&](int y) {
	//for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			// 乱数
			Random rng(y * width + x);

			// ピクセル色の初期化
			Color &pixel = image[y * width + x];
			pixel = Color(0.0, 0.0, 0.0);

			// サンプリング
			for (int s = 0; s < samples; s++) {
				// テントフィルターによってサンプリング
				const double r1 = 2.0 * rng.next01();
				const double r2 = 2.0 * rng.next01();
				const double dx = r1 < 1.0 ? sqrt(r1) - 1.0 : 1.0 - sqrt(2.0 - r1);
				const double dy = r2 < 1.0 ? sqrt(r2) - 1.0 : 1.0 - sqrt(2.0 - r2);

				const double px = (x + dx + 0.5) / width - 0.5;
				const double py = ((height - y - 1) + dy + 0.5) / height - 0.5;

				// 放射輝度の計算
				const Vec dir = cx * px + cy * py + camera.dir;
				const Ray ray(camera.org + dir * 13.0, Normalize(dir));
				const Color L = radiance(ray, Medium(), rng, 0, maxDepth);
				Assertion(L.isValid(), "Radiance is invalid: (%f, %f %f)", L.x, L.y, L.z);

				pixel = pixel + L;
			}
			pixel = pixel / samples;

			// 進行度の表示
			mtx.lock();
			progressBar(++progress, width * height);
			mtx.unlock();
			}
		});
	clock_t end = clock();
	std::cout << (float)(end - start) / CLOCKS_PER_SEC << std::endl;
	// PPMファイルを保存
	const std::string outfile = std::string(OUTPUT_DIRECTORY) + "image.ppm";
	save_ppm_file(outfile, image.get(), width, height);
}
