#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector>

namespace rt {

const double PI = 3.14159265358979323846;
const double INF = 1e20;
const double EPS = 1e-6;
const double MaxDepth = 5;

struct Vec {
	double x, y, z;
	Vec(const double x_ = 0, const double y_ = 0, const double z_ = 0) : x(x_), y(y_), z(z_) {}
	inline Vec operator+(const Vec &b) const {return Vec(x + b.x, y + b.y, z + b.z);}
	inline Vec operator-(const Vec &b) const {return Vec(x - b.x, y - b.y, z - b.z);}
	inline Vec operator*(const double b) const {return Vec(x * b, y * b, z * b);}
	inline Vec operator/(const double b) const {return Vec(x / b, y / b, z / b);}
	inline const double LengthSquared() const { return x*x + y*y + z*z; }
	inline const double Length() const { return sqrt(LengthSquared()); }
};
inline Vec operator*(double f, const Vec &v) { return v * f; }
inline Vec Normalize(const Vec &v) { return v / v.Length(); }
inline const Vec Multiply(const Vec &v1, const Vec &v2) {
	return Vec(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}
inline const double Dot(const Vec &v1, const Vec &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
inline const Vec Cross(const Vec &v1, const Vec &v2) {
	return Vec((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
}
typedef Vec Color;

struct Ray {
	Vec org, dir;
	Ray(const Vec org_, const Vec &dir_) : org(org_), dir(dir_) {}
};

struct Material {
	double diffuse;
	Color color;
	double highlight;
	bool mirror;
	int texture;

	Material(double diffuse_, Color &color_, double highlight_, bool mirror_, int texture_) :
	diffuse(diffuse_), color(color_), highlight(highlight_), mirror(mirror_), texture(texture_) {
	}
};

class Object {
protected:
	double reverse;

public:
	Vec position;
	Material mat;
	Object(const Material &mat_, const Vec &position_, const double reverse_) :
	mat(mat_), position(position_), reverse(reverse_) {
	}

	virtual bool in(const Vec &p) = 0;
	virtual bool intersect(const Ray &ray, double *near, double *far, Vec *near_normal, Vec *far_normal) = 0;
};

inline double sign(const double x) {
	return x > 0.0 ? 1.0 : -1.0;
}

struct Quadric : public Object {
private:
	Vec param;
	double A, B, C;
	double A2, B2, C2;
	double signA, signB, signC;

	double a, b, c;
public:
	Quadric(const Material &mat_, const Vec &position_,  const double reverse_, const Vec &param_) : 
	  Object(mat_, position_, reverse_), param(param_) {
		  A = param.x; B = param.y; C = param.z;
		  A2 = A * A; B2 = B * B; C2 = C * C;
		  signA = sign(A); signB = sign(B); signC = sign(C);
		  a = A2 == 0 ? 0.0 : signA / A2;
		  b = B2 == 0 ? 0.0 : signB / B2; 
		  c = C2 == 0 ? 0.0 : signC / C2;
	  }
 
	bool intersect(const Ray &ray, double *near, double *far, Vec *near_normal, Vec *far_normal) {
		*near = INF;
		*far  = 0.0;

		const Vec o = ray.org - position;

		const double X = a * ray.dir.x * ray.dir.x + b * ray.dir.y * ray.dir.y + c * ray.dir.z * ray.dir.z;
		const double Y = 2.0 * a * o.x * ray.dir.x + 2.0 * b * o.y * ray.dir.y + 2.0 * c * o.z * ray.dir.z;
		const double Z = a * o.x * o.x + b * o.y * o.y + c * o.z * o.z - 1.0;
		
		if (fabs(X) <= EPS)
			return false;
		
		const double W = Y * Y - 4.0 * X * Z;
		if (W < 0.0)
			return false;

		const double t0 = (-Y - sqrt(W)) / (2.0 * X); // near
		const double t1 = (-Y + sqrt(W)) / (2.0 * X); // far

		const Vec p0 = o + t0 * ray.dir;
		const Vec p1 = o + t1 * ray.dir;

		if (t1 < 0.0)
			return false;
		if (t0 < 0.0) {
			*near = t1;
			*near_normal = reverse * Normalize(Vec(2.0 * a * p1.x, 2.0 * b * p1.y, 2.0 * c * p1.z));
			return true;
		}
		*near = t0;
		*far  = t1;
		*near_normal = reverse * Normalize(Vec(2.0 * a * p0.x, 2.0 * b * p0.y, 2.0 * c * p0.z));
		*far_normal  = reverse * Normalize(Vec(2.0 * a * p1.x, 2.0 * b * p1.y, 2.0 * c * p1.z));

		return true;
			
	}

	bool in_sub(const Vec &p) {
		const double d = a * p.x * p.x + b * p.y * p.y + c * p.z * p.z - 1.0;

		return (reverse * d) < 0.0;
	}
	bool in(const Vec &p) {
		bool ret = in_sub(p - position);
		return ret;
	}
};

struct Plane : public Object {
private:
	Vec normal;
public:
	Plane(const Material &mat_, const Vec &position_,  const double reverse_, const Vec &normal_) : 
	  Object(mat_, position_, reverse_), normal(normal_) {
	  }
	  
	bool intersect(const Ray &ray, double *near, double *far, Vec *near_normal, Vec *far_normal) {
		*near = INF;
		*far  = 0.0;

		const double denom = Dot(ray.dir, normal);
		if (fabs(denom) <= EPS) {
			return false;
		} else {
			const double t = -1.0 * Dot(ray.org - position, normal) / denom;
			if (t < 0.0)
				return false;
			*near = t;
			*near_normal = reverse * normal;
			return true;
		}
	}

	bool in_sub(const Vec &p) {
		const Vec d = p - position;
		return reverse * Dot(d, normal) < 0.0;
	}
	bool in(const Vec &p) {
		bool ret = in_sub(p);
		return ret;
	}
};

struct Box : public Object {
private:
	Vec size;
public:
	Box(const Material &mat_, const Vec &position_,  const double reverse_, const Vec &size_) : 
	  Object(mat_, position_, reverse_), size(size_) {

		  position = position - size;
		  size = 2.0 * size;
	  }
	  
	bool intersect(const Ray &ray, double *near, double *far, Vec *near_normal, Vec *far_normal) {
		*near = INF;
		*far  = 0.0;
		// z = position.z
		const double t0 = (position.z - ray.org.z) / ray.dir.z;
		const double x0 = ray.org.x + t0 * ray.dir.x;
		const double y0 = ray.org.y + t0 * ray.dir.y;
		if (position.x <= x0 && x0 <= position.x + size.x && 
			position.y <= y0 && y0 <= position.y + size.y && t0 > 0.0) {
				if (t0 < *near) {
					*near = t0;
					*near_normal = Vec(0.0, 0.0, -1.0);
				}
				if (*far < t0) {
					*far  = t0;
					*far_normal = Vec(0.0, 0.0, -1.0);
				}
		}
		// z = positn.z + size.z
		const double t1 = (position.z + size.z - ray.org.z) / ray.dir.z;
		const double x1 = ray.org.x + t1 * ray.dir.x;
		const double y1 = ray.org.y + t1 * ray.dir.y;
		if (position.x <= x1 && x1 <= position.x + size.x && 
			position.y <= y1 && y1 <= position.y + size.y && t1 > 0.0) {
				if (t1 < *near) {
					*near = t1;
					*near_normal = Vec(0.0, 0.0, 1.0);
				}
				if (*far < t1) {
					*far  = t1;
					*far_normal = Vec(0.0, 0.0, 1.0);
				}
		}
		// x = position.x
		const double t2 = (position.x - ray.org.x) / ray.dir.x;
		const double y2 = ray.org.y + t2 * ray.dir.y;
		const double z2 = ray.org.z + t2 * ray.dir.z;
		if (position.y <= y2 && y2 <= position.y + size.y && 
			position.z <= z2 && z2 <= position.z + size.z && t2 > 0.0) {
				if (t2 < *near) {
					*near = t2;
					*near_normal = Vec(-1.0, 0.0, 0.0);
				}
				if (*far < t2) {
					*far  = t2;
					*far_normal = Vec(-1.0, 0.0, 0.0);
				}
		}
		// x = position.x + size.x
		const double t3 = (position.x + size.x - ray.org.x) / ray.dir.x;
		const double y3 = ray.org.y + t3 * ray.dir.y;
		const double z3 = ray.org.z + t3 * ray.dir.z;
		if (position.y <= y3 && y3 <= position.y + size.y && 
			position.z <= z3 && z3 <= position.z + size.z && t3 > 0.0) {
				if (t3 < *near) {
					*near = t3;
					*near_normal = Vec(1.0, 0.0, 0.0);
				}
				if (*far < t3) {
					*far  = t3;
					*far_normal = Vec(1.0, 0.0, 0.0);
				}
		}
		// y = position.y
		const double t4 = (position.y - ray.org.y) / ray.dir.y;
		const double x4 = ray.org.x + t4 * ray.dir.x;
		const double z4 = ray.org.z + t4 * ray.dir.z;
		if (position.x <= x4 && x4 <= position.x + size.x && 
			position.z <= z4 && z4 <= position.z + size.z && t4 > 0.0) {
				if (t4 < *near) {
					*near = t4;
					*near_normal = Vec(0.0, -1.0, 0.0);
				}
				if (*far < t4) {
					*far  = t4;
					*far_normal = Vec(0.0, -1.0, 0.0);
				}
		}
		// x = position.x + size.x
		const double t5 = (position.y + size.y - ray.org.y) / ray.dir.y;
		const double x5 = ray.org.x + t5 * ray.dir.x;
		const double z5 = ray.org.z + t5 * ray.dir.z;
		if (position.x <= x5 && x5 <= position.x + size.x && 
			position.z <= z5 && z5 <= position.z + size.z && t5 > 0.0) {
				if (t5 < *near) {
					*near = t5;
					*near_normal = Vec(0.0, 1.0, 0.0);
				}
				if (*far < t5) {
					*far  = t5;
					*far_normal = Vec(0.0, 1.0, 0.0);
				}
		}

		if (*near >= INF && 0 >= *far)
			return false;
		else {
			*near_normal = reverse * (*near_normal);
			*far_normal  = reverse * (*far_normal);
			return true;
		}
	}

	bool in_sub(const Vec &p) {
		if (p.x < position.x || position.x + size.x < p.x)
			return false;
		if (p.y < position.y || position.y + size.y < p.y)
			return false;
		if (p.z < position.z || position.z + size.z < p.z)
			return false;

		return true;
	}
	bool in(const Vec &p) {
		bool ret = in_sub(p);
		if (reverse < 0.0)
			return !ret;
		return ret;
	}
};


// *** .hdrフォーマットで出力するための関数 ***
struct HDRPixel {
	unsigned char r, g, b, e;
	HDRPixel(const unsigned char r_ = 0, const unsigned char g_ = 0, const unsigned char b_ = 0, const unsigned char e_ = 0) :
	r(r_), g(g_), b(b_), e(e_) {};
	unsigned char get(int idx) {
		switch (idx) {
		case 0: return r;
		case 1: return g;
		case 2: return b;
		case 3: return e;
		} return 0;
	}

};

// doubleのRGB要素を.hdrフォーマット用に変換
HDRPixel get_hdr_pixel(const Color &color) {
	double d = std::max(color.x, std::max(color.y, color.z));
	if (d <= 1e-32)
		return HDRPixel();
	int e;
	double m = frexp(d, &e); // d = m * 2^e
	d = m * 256.0 / d;
	return HDRPixel(color.x * d, color.y * d, color.z * d, e + 128);
}

// 書き出し用関数
void save_hdr_file(const std::string &filename, const Color* image, const int width, const int height) {
	FILE *fp = fopen(filename.c_str(), "wb");
	if (fp == NULL) {
		std::cerr << "Error: " << filename << std::endl;
		return;
	}
	// .hdrフォーマットに従ってデータを書きだす
	// ヘッダ
	unsigned char ret = 0x0a;
	fprintf(fp, "#?RADIANCE%c", (unsigned char)ret);
	fprintf(fp, "# Made with 100%% pure HDR Shop%c", ret);
	fprintf(fp, "FORMAT=32-bit_rle_rgbe%c", ret);
	fprintf(fp, "EXPOSURE=1.0000000000000%c%c", ret, ret);

	// 輝度値書き出し
	fprintf(fp, "-Y %d +X %d%c", height, width, ret);
	for (int i = height - 1; i >= 0; i --) {
		std::vector<HDRPixel> line;
		for (int j = 0; j < width; j ++) {
			HDRPixel p = get_hdr_pixel(image[j + i * width]);
			line.push_back(p);
		}
		fprintf(fp, "%c%c", 0x02, 0x02);
		fprintf(fp, "%c%c", (width >> 8) & 0xFF, width & 0xFF);
		for (int i = 0; i < 4; i ++) {
			for (int cursor = 0; cursor < width;) {
				const int cursor_move = std::min(127, width - cursor);
				fprintf(fp, "%c", cursor_move);
				for (int j = cursor; j < cursor + cursor_move; j ++)
					fprintf(fp, "%c", line[j].get(i));
				cursor += cursor_move;
			}
		}
	}

	fclose(fp);
}

class Scene {
public:
	Vec screen_center;
	double eye0, eye1;
	int light_num;
	double light0, light1;
	double light_intensity;

	std::vector<Object*> objects;

	struct AndNetwork {
		std::vector<int> object_id;
	};

	struct OrNetwork {
		int range_primitive_id;
		std::vector<int> and_network_id;
	};

	std::vector<AndNetwork> and_networks;
	std::vector<OrNetwork> or_networks;

	void load(const char* filename) {
		std::vector<float> scene;
		FILE *fp = fopen(filename, "rt");
		if (fp != NULL) {
			for (;;) {
				float v;
				int ret = fscanf(fp, "%f", &v);
				if (ret == -1)
					break;
				scene.push_back(v);
			}
			fclose(fp);
		}

		int index = 0;

		screen_center.x = scene[index++];
		screen_center.y = scene[index++];
		screen_center.z = scene[index++];

		eye0 = scene[index++] * PI / 180.0;
		eye1 = scene[index++] * PI / 180.0;

		light_num = scene[index++];
		light0 = scene[index++] * PI / 180.0;
		light1 = scene[index++] * PI / 180.0;
		light_intensity = scene[index++];

		for (;;) {
			if (scene[index] == -1.0)
				break;

			const int texture = scene[index++];
			const int shape = scene[index++];
			const bool mirror = scene[index++] == 2 ? true : false;
			const int rotation = scene[index++];
			const double sx = scene[index++];
			const double sy = scene[index++];
			const double sz = scene[index++];
			const double px = scene[index++];
			const double py = scene[index++];
			const double pz = scene[index++];
			const double reverse = scene[index++];
			const double diffuse = scene[index++];
			const double highlight = scene[index++];
			const double R = scene[index++];
			const double G = scene[index++];
			const double B = scene[index++];

			const Material mat(diffuse, Color(R, G, B), highlight == 2, mirror, texture);
			switch (shape) {
			case 1: 
				{
					objects.push_back(new Box(mat, Vec(px, py, pz), reverse, Vec(sx, sy, sz)));
				}
				break;
			case 2: 
				{
					objects.push_back(new Plane(mat, Vec(px, py, pz), reverse, Vec(sx, sy, sz)));
				}
				break;
			case 3: 
				{
					objects.push_back(new Quadric(mat, Vec(px, py, pz), reverse, Vec(sx, sy, sz)));
				}
				break;
			}
		}
		index ++;

		// ANDネットワーク
		for (;;) {
			if (scene[index] == -1.0)
				break;

			AndNetwork n;
			for (;;) {
				if (scene[index] == -1.0)
					break;
				int a = scene[index ++];
				n.object_id.push_back(a);
			}
			and_networks.push_back(n);
			index ++;
		}
		index ++;
		
		// ORネットワーク
		for (;;) {
			if (scene[index] == -1.0)
				break;

			OrNetwork n;
			int range_primitive_id = scene[index++];
			n.range_primitive_id = range_primitive_id;
			for (;;) {
				if (scene[index] == -1.0)
					break;
				int a = scene[index ++];
				n.and_network_id.push_back(a);
			}
			or_networks.push_back(n);
			index ++;
		}
		index ++;
	}

	bool intersect(const Ray &ray, double *t, Vec *normal, int *obj_id) {
		double t_ = INF;
		int obj_id_;
		Vec normal_;
		for (int i = 0; i < or_networks.size(); i++) {
			OrNetwork& or = or_networks[i];

			double a_, b_;
			Vec c_, d_;
			if (or.range_primitive_id != 99 && !objects[or.range_primitive_id]->intersect(ray, &a_, &b_, &c_, &d_))
				continue;

			for (int j = 0; j < or.and_network_id.size(); j ++) {
				AndNetwork& and = and_networks[or.and_network_id[j]];

				for (int k = 0; k < and.object_id.size(); k ++) {
					double near = INF, far = 0.0;
					Vec nearn, farn;
					if (objects[and.object_id[k]]->intersect(ray, &near, &far, &nearn, &farn)) {
						bool flag = true;
						const Vec p = ray.org + near * ray.dir;
						// 交差点がほかのオブジェクトの内部にあるかどうかを判定
						for (int n = 0; n < and.object_id.size(); n ++) {
							if (k == n)
								continue;
							if (!objects[and.object_id[n]]->in(p)) {
								flag = false;
								break;
							}
						}
						// OK
						if (flag) {
							if (near < t_) {
								t_ = near;
								obj_id_ = and.object_id[k];
								normal_ = nearn;
							}
						}

						if (far > 0.0) {
							bool flag = true;
							const Vec p = ray.org + far * ray.dir;
							// 交差点がほかのオブジェクトの内部にあるかどうかを判定
							for (int n = 0; n < and.object_id.size(); n ++) {
								if (k == n)
									continue;
								if (!objects[and.object_id[n]]->in(p)) {
									flag = false;
									break;
								}
							}
							// OK
							if (flag) {
								if (far < t_) {
									t_ = far;
									obj_id_ = and.object_id[k];
									normal_ = farn;
								}
							}
						}
					}
				}
			}
		}

		if (t_ < INF) {
			*t = t_;
			*normal = normal_;
			*obj_id = obj_id_;
			return true;
		}

		return false;
	}
	inline double rand01() { return (double)rand()/RAND_MAX; }

	Color radiance(const Ray &ray, const int depth) {
		double t = INF;
		int obj_id = -1;
		Vec normal;

		if (!intersect(ray, &t, &normal, &obj_id))
			return Color();
		const Vec orienting_normal = Dot(normal, ray.dir) < 0.0 ? normal : -1.0 * normal;

		const Vec hitpoint = ray.org + (t - EPS) * ray.dir;
		double russian_roulette_probability = 0.5;

		if (depth > MaxDepth) {
			if (rand01() >= russian_roulette_probability)
				return Color();
		} else
			russian_roulette_probability = 1.0; // ロシアンルーレット実行しなかった
		
		if (objects[obj_id]->mat.mirror) {
			return radiance(Ray(hitpoint, ray.dir - normal * 2.0 * Dot(normal, ray.dir)), depth + 1)
				/ russian_roulette_probability;
		}

		if (objects[obj_id]->mat.texture == 9) {
			return objects[obj_id]->mat.color / 255.0;
		}

		// orienting_orienting_normalの方向を基準とした正規直交基底(w, u, v)を作る。この基底に対する半球内で次のレイを飛ばす
		Vec w, u, v;
		w = orienting_normal;
		if (fabs(w.x) > 0.1)
			u = Normalize(Cross(Vec(0.0, 1.0, 0.0), w));
		else
			u = Normalize(Cross(Vec(1.0, 0.0, 0.0), w));
		v = Cross(w, u);
		// コサイン項を使った重点的サンプリング
		const double r1 = 2 * PI * rand01();
		const double r2 = rand01(), r2s = sqrt(r2);
		Vec dir = Normalize((u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1.0 - r2)));

		// チェッカー
		Color diffuse_color = objects[obj_id]->mat.color;
		if (objects[obj_id]->mat.texture == 1) {
			const Vec diff = hitpoint - objects[obj_id]->position;

			const double w0 = floor(diff.x * 0.05) * 20.0;
			const double w1 = floor(diff.z * 0.05) * 20.0;
			const int flag0 = (diff.x - w0) < 10.0;
			const int flag1 = (diff.z - w1) < 10.0;
			if (flag0) {
				if (flag1)  {
					diffuse_color.y = 255;
				} else {
					diffuse_color.y = 0;
				}
			} else {
				if (flag1)  {
					diffuse_color.y = 0;
				} else {
					diffuse_color.y = 255;
				}
			}
			
			return Multiply(objects[obj_id]->mat.diffuse * diffuse_color / 255.0, radiance(Ray(hitpoint, dir), depth + 1))
				/ russian_roulette_probability;
		} else if (objects[obj_id]->mat.texture == 2) {
			const double w2 = pow(sin(hitpoint.y * 0.25), 2);
			diffuse_color.x = 255.0 * w2;
			diffuse_color.y = 255.0 * (1.0 - w2);
			
			return Multiply(objects[obj_id]->mat.diffuse * diffuse_color / 255.0, radiance(Ray(hitpoint, dir), depth + 1))
				/ russian_roulette_probability;
		} 
		// 各種パラメータ。
		// 決め方がださいね。

		double ro = 1.0;
		double alpha_x = 0, alpha_y = 0;

		alpha_x = 0.4;
		alpha_y = 0.4;
		const Vec in = -1.0 * ray.dir;
		Vec halfv;

		// 重点サンプリングする。ただし、生成されるハーフベクトルしだいでは
		// 反射方向が半球外に出てしまうのでそういう場合は棄却する。
		do {
			const double u1 = rand01();
			const double u2 = rand01();
			double phi = atan(alpha_y / alpha_x * tan(2.0 * PI * u2));
			if (0.25 <= u2 && u2 <= 0.75)
				phi += PI;
			else if (0.75 < u2)
				phi += 2.0 * PI;
			const double theta = atan(sqrt(-log(u1) / (pow(cos(phi), 2) / pow(alpha_x, 2) + pow(sin(phi), 2) / pow(alpha_y, 2))));

			halfv = Normalize((u * cos(phi) * sin(theta) + v * sin(phi) * sin(theta) + w * cos(theta)));
			dir = 2.0 * Dot(in, halfv) * halfv - in;
		} while (Dot(orienting_normal, dir) < 0.0);

		// 重み。brdf * cosθ / pdf(ω)をまとめたものになっている。
		const double weight = ro * Dot(halfv, in) * pow(Dot(halfv, orienting_normal), 3) *
			sqrt(Dot(dir, orienting_normal) / Dot(in, orienting_normal));


		return Multiply(objects[obj_id]->mat.diffuse * diffuse_color / 255.0, radiance(Ray(hitpoint, dir), depth+1))
			* weight
			/ russian_roulette_probability;

	}

	
	// concentricにサンプリング
	void concentric_sample_disk(double u1, double u2, double *dx, double *dy) {
		double r, theta;
		// [0, 1]の一様乱数u1,u2を[-1, 1]の一様乱数sx,syに写像
		const double sx = 2 * u1 - 1;
		const double sy = 2 * u2 - 1;


		// sx, syが0,0だった場合は特別に処理
		if (sx == 0.0 && sy == 0.0) {
			*dx = 0.0;
			*dy = 0.0;
			return;
		}
	// 四つに分割した円の各部位で別々の処理になる
		if (sx >= -sy) {
			if (sx > sy) {
				r = sx;
				if (sy > 0.0) theta = sy/r;
				else theta = 8.0f + sy/r;
			}
			else {
				r = sy;
				theta = 2.0f - sx/r;
			}
		}
		else {
			if (sx <= sy) {
				r = -sx;
				theta = 4.0f - sy/r;
			}
			else {
				r = -sy;
				theta = 6.0f + sx/r;
			}
		}
		theta *= PI / 4.f;
		*dx = r * cosf(theta);
		*dy = r * sinf(theta);
	}

	void render(const int width, const int height) {
		const int image_center_x = width / 2;
		const int image_center_y = height / 2;
		const double scan_pitch = (double)128.0 / width;
		const Vec screen_z_dir(cos(eye0) * sin(eye1) * 200.0, sin(eye0) * -200.0, cos(eye0) * cos(eye1) * 200.0);
		const Vec screen_x_dir(cos(eye1), 0.0, -sin(eye1));
		const Vec screen_y_dir(-sin(eye0)*sin(eye1), -cos(eye0), -sin(eye0)*cos(eye1));
		const Vec viewpoint = screen_center - screen_z_dir;

		Color* image = new Color[width * height];

		const int samples = 256;

//#pragma omp parallel for schedule(dynamic, 1) num_threads(3)
		for (int iy = 0; iy < height; iy ++) {
			srand(iy * iy + iy);
			std::cout << iy << std::endl;
			for (int ix = 0; ix < width; ix ++) {
				Color accum;
				for (int sy = 0; sy < 2; sy ++) {
					for (int sx = 0; sx < 2; sx ++) {
						const double rx = sx / 2.0;
						const double ry = sy / 2.0;
						const double ydisp = scan_pitch * ((iy + ry) - image_center_y);
						const double xdisp = scan_pitch * ((ix + rx) - image_center_x);
						const Vec dir = xdisp * screen_x_dir + ydisp * screen_y_dir + screen_z_dir;
						Ray ray(viewpoint, Normalize(dir));
						
						for (int s = 0; s < samples; s ++) {
							accum = accum + radiance(ray, 0);
						}

					}
				}
				image[(height - iy - 1) * width + ix] = accum / samples / 4;
			}
		}

		save_hdr_file("hoge.hdr", image, width, height);

		delete[] image;
	}
};


};

int main() {
	rt::Scene s;
	s.load("test.sld");

	s.render(128, 128);

	return 0;
}